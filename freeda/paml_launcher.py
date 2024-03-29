#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copyright 2022 - Damian Dudka and R. Brian Akins - contact: damiandudka0@gmail.com

This file is part of FREEDA.

FREEDA is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

FREEDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with FREEDA.
If not, see <https://www.gnu.org/licenses/>.

"""


"""

Analyses the final cds, gets a gene tree based on translated cds and runs PAML.

"""

from freeda import fasta_reader
from freeda import genomes_preprocessing
from freeda import control_file
from freeda import gui_logging_handler
from freeda import pyinstaller_compatibility
from Bio.Align.Applications import PrankCommandline
from Bio import pairwise2
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML.chi2 import cdf_chi2

import datetime
import time
import glob
import os
import re
import logging
import shutil
import subprocess
import copy
import pyensembl


def analyze_final_cds(wdir, ref_species, result_path, all_genes, codon_frequencies, gui=None, logging_window=None,
                      final_excluded_species=None, subgroup=None):
    """Main function controlling building final alignments, gene trees and PAML anaylsis"""

    failed_paml = []
    start_time = time.time()
    day = datetime.datetime.now().strftime("-%m-%d-%Y-%H-%M")
    PAML_logfile_name = "PAML" + day + ".log"

    # initiate log file to record PAML analysis by reseting the handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    if gui:
        # make a new handler
        text_handler = gui_logging_handler.TextHandler(logging_window)
        # configure the logger
        logging.basicConfig(filename=PAML_logfile_name, level=logging.INFO, format="%(message)s")
        logger = logging.getLogger()
        logger.addHandler(text_handler)

    else:
        # configure the logger
        logging.basicConfig(filename=PAML_logfile_name, level=logging.INFO, format="%(message)s")

    logging.info("\n\n======================== PAML ANALYSIS ========================")

    # make an empty template for all possible species
    all_species = [names[0] for names in genomes_preprocessing.get_names(wdir, ref_species, final_excluded_species, subgroup)]

    all_species_dict = {ref_species: ""}
    for species in all_species:
        all_species_dict[species] = ""
    
    nr_of_species_total_dict = {}
    genes_under_pos_sel = {"F3X4": [], "F61": []}

    for gene in all_genes:
        
        # check if this gene was already analyzed
        if os.path.isdir(result_path + gene + "/" + "PAML_" + gene + "_" + "F3X4") \
                or os.path.isdir(result_path + gene + "/" + "PAML_" + gene + "_" + "F61"):
            message = "\n################\n\n PAML analysis has been already performed for : %s (skipping)" % gene
            print(message)
            logging.info(message)
           
            # get how many species were analyzed
            path_to_final_species = result_path + gene + "/" + gene + "_final.fasta"
            with open(path_to_final_species, "r") as f:
                nr_of_species = f.read().count(">")
                nr_of_species_total_dict[gene] = nr_of_species
           
        # otherwise proceed with the analysis
        else:
            # need to use deepcopy function to make an actual dictionary copy
            final_species = copy.deepcopy(all_species_dict)
            final_species_headers = []
        
            # read each cds fasta file and put them into the empty final_species dict
            cds_file_path = result_path + gene + ".fasta"
        
            if os.path.isfile(cds_file_path) is False:
                message = "\nThis file has not been found : %s." % cds_file_path
                print(message)
                logging.info(message)
                continue
        
            with open(cds_file_path, "r") as f:
                file = f.readlines()
            
                # count lines
                line_nr = 0
                for line in file:
                    line_nr += 1
                
                    # first line is always the ref species header
                    if line_nr == 1:
                        ref_head = line.split(gene + "_")[1].rstrip("\n")
                        continue
                
                    # second line is always the ref species sequence
                    if line_nr == 2:
                        ref_seq = line.rstrip("\n")
                        final_species[ref_head] = ref_seq
                        continue
                
                    # next lines are headers and sequences of cloned cds
                    if line.startswith(">"):
                        head = line.split(gene + "_")[1].rstrip("\n")
                        continue
                
                    else:
                        seq = line.rstrip("\n")
                        seq_no_dashes = seq.replace("-", "")
                        if len(seq_no_dashes) / len(ref_seq) >= 0.90:
                            final_species[head] = seq
                            final_species_headers.append(head)
               
            # log species
            final_species_headers.insert(0, ref_species)
            nr_of_species = len(final_species_headers)

            message = "\n\n --------- * %s * --------- \n\n" % gene
            print(message)
            logging.info(message)
        
            # write all sequences into a file by reading through the dictionary
            final_cds_file = gene + "_final.fasta"
            with open(final_cds_file, "w") as f:
                for species in final_species:
                    
                    if final_species[species] == "":
                        continue
                    else:
                        f.write(">" + species + "\n")
                        f.write(final_species[species] + "\n")

            gene_folder_path = result_path + gene
            
            shutil.move(cds_file_path, gene_folder_path)
            shutil.move(final_cds_file, gene_folder_path)

            # make and copy control_file from working directory into gene_folder_path
            if not os.path.isfile(wdir + "control_file_F3X4.ctl") \
                    or not os.path.isfile(wdir + "control_file_F61.ctl"):
                control_file.make_control_file(wdir)

            # define which codon frequency is used - allow only F3x4 (default) or F61
            if str(codon_frequencies).upper().replace(" ", "") != "F61" \
                and str(codon_frequencies).upper().replace(" ", "") != "F3X4" \
                and str(codon_frequencies).upper().replace(" ", "") != "F3X4,F61" \
                and str(codon_frequencies).upper().replace(" ", "") != "F61,F3X4":
                codon_frequencies = "F3X4"

            # align the final cds sequences
            out_msa, aligner = align_final_cds(wdir, gene, result_path)

            # write the alignment into a file from dict and swap full names
            path = gene_folder_path + "/" + out_msa
            all_seq_dict = fasta_reader.alignment_file_to_dict(wdir, ref_species, path, raw_cds=True)
            new_result_path = result_path.replace("Raw_data/", "Results/") + "Nucleotide_alignments/" + out_msa
            with open(new_result_path, 'w') as file:
                for ab, s in all_seq_dict.items():
                    file.write(ab + "\n")
                    file.write(s + "\n")
            genomes_preprocessing.substitute_abbreviations(ref_species, new_result_path)

            # check and eliminate insertions that cause dashes in ref species
            correction, corrected_filename = eliminate_all_insertions(gene_folder_path, out_msa)
            if correction is True:
                side_note = " \n---- Insertions detected in " + gene + " alignment -> " \
                    "these bp positions were removed in all species forcing conserved alignment\n"
                print(side_note)
                logging.info(side_note)
                no_dashes_out_msa = corrected_filename
            else:
                no_dashes_out_msa = out_msa

            # remove STOP codons
            final_cds_file_no_STOP = STOP_remover(gene_folder_path, no_dashes_out_msa, gene, aligner)

            # check for rare frameshifts in cloned cds
            to_delete = cloned_cds_frameshift_checkpoint(wdir, ref_species, gene, final_cds_file_no_STOP)
            shutil.move(final_cds_file_no_STOP, gene_folder_path)

            # run gBLOCK
            raw_out_Gblocks_filename = run_Gblocks(final_cds_file_no_STOP, gene, result_path, aligner)
            # Gblocks will fail if not enough species in alignment -> go to next gene
            if raw_out_Gblocks_filename is None:
                failed_paml.append(gene)
                continue

            shutil.move(wdir + raw_out_Gblocks_filename, gene_folder_path)

            # double check for artificial STOP codons introduced by Gblocks -> force conserved alignment
            # it overwrites the previous out_Gblocks file (same name)
            post_Gblocks_STOP_remover(gene_folder_path, raw_out_Gblocks_filename, gene)

            # translate all sequences (inside of the gene folder)
            raw_translated_path, filepath_to_translate = translate_Gblocks(wdir, gene_folder_path,
                                                                           raw_out_Gblocks_filename, gene, ref_species)

            # check for frameshifts in gene seq and eliminate species from gene and cds
            translated_path, out_Gblocks, to_delete = eliminate_frameshits_cds(wdir, ref_species, gene,
                                                                               raw_translated_path,
                                                                               filepath_to_translate,
                                                                               to_delete, aligner)

            # check compatibility of most distant cloned seq with its known gene sequence (e.g. rat for mouse)
            check_compatibility(ref_species, gene, translated_path)

            # correct number of species in the final analysis
            final_nr_of_species = nr_of_species - len(to_delete)
            nr_of_species_total_dict[gene] = final_nr_of_species
            final_species_headers_corrected = [species for species in final_species_headers if ">" + species not in to_delete]
            message = "\n Final species cloned and aligned (+ ref) for %s : %s %s \n" \
                      % (gene, str(final_nr_of_species), str(final_species_headers_corrected))
            print(message)
            logging.info(message)

            # convert fasta alignment to phylip format
            convert_to_phylip(out_Gblocks)

            # run RAxML (and move all the RAxML files to gene folder)
            try:
                best_tree_path = run_RAxML(gene, gene_folder_path, out_Gblocks)

            except FileNotFoundError:
                failed_paml.append(gene)
                message = "...PROBLEM... : The best gene tree was NOT found for %s " \
                          "-> likely too few cds pass coverage threshold (90 percent default)" % gene
                print(message)
                logging.info(message)
                continue

            for codon_frequency in codon_frequencies.split(","):
                codon_frequency = codon_frequency.upper().replace(" ", "")  # get rid of the spaces if present

                # generate a PAML folder for a given gene
                PAML_path = gene_folder_path + "/PAML_" + gene + "_" + codon_frequency
                os.makedirs(PAML_path)

                control_file_name = "control_file_" + codon_frequency + ".ctl"

                # copy all required files into PAML path
                shutil.copy(control_file_name, PAML_path + "/" + control_file_name)
                shutil.copy(out_Gblocks, PAML_path + "/input.phy")
                shutil.copy(best_tree_path, PAML_path + "/gene.tree")

                # rename path
                path = result_path.replace("Raw_data/", "Results/Gene_trees/") + gene + ".tree"
                shutil.copy(best_tree_path, path)

                # match abbreviations to species names
                genomes_preprocessing.substitute_abbreviations(ref_species, path, tree=True)

                # rename and copy the protein alignment into results (not final at this point yet)
                shutil.copy(translated_path, result_path + gene + "_protein_alignment.fasta")

                # run PAML
                message = "\n.........Running PAML (%s) for gene: %s.........\n" % (codon_frequency, gene)
                print(message)
                logging.info(message)

                M2a_M1a, M8_M7 = run_PAML(wdir, gene, PAML_path, control_file_name, codon_frequency)

                if M8_M7 < 0.05:
                    genes_under_pos_sel[codon_frequency].append(gene)

                message = "\n -> PAML p-values (%s) for gene %s : M2a v M1a - %s and M8 v M7 - %s\n" \
                    % (codon_frequency, gene, str(M2a_M1a), str(M8_M7))
                print(message)
                logging.info(message)

    shutil.move(wdir + PAML_logfile_name, result_path)

    # mark the end of the analysis
    message = ("\n --------------->  PAML analysis completed in %s minutes or %s hours"
               % ((time.time() - start_time)/60, (time.time() - start_time)/60/60))
    print(message)
    logging.info(message)

    return nr_of_species_total_dict, PAML_logfile_name, day, failed_paml, genes_under_pos_sel


def check_compatibility(ref_species, gene, translated_path):
    """Checks cloned protein sequence from opposite species against ensembl entry (e.g. rat for mouse)"""

    if ref_species == "Mm":
        species = "Rn"
        # get known protein sequence from the most distant species
        distant_seq = get_most_distant_prot_seq(species, gene)

    elif ref_species == "Rn":
        species = "Mm"
        # get known protein sequence from the most distant species
        distant_seq = get_most_distant_prot_seq(species, gene)

    elif ref_species == "Fc":
        species = "Cf"
        # get known protein sequence from the most distant species
        distant_seq = get_most_distant_prot_seq(species, gene)

    elif ref_species == "Cf":
        species = "Fc"
        # get known protein sequence from the most distant species
        distant_seq = get_most_distant_prot_seq(species, gene)

    # NOT SUPPORTING PRIMATES AND BIRDS
    else:
        return

    # get cloned protein sequence from the most distant species
    record_sequence = False
    seq = ""
    with open(translated_path, "r") as f:
        file = f.readlines()
        for line in file:
            if line.startswith(">" + species):
                record_sequence = True
                continue
            if record_sequence:
                seq = line.rstrip("\n")
                break

    # align them and write into dicts
    aln = pairwise2.align.globalxs(distant_seq, seq, open=-0.5, extend=-0.1)

    if not distant_seq and not seq:
        message = "\n...WARNING... : Gene %s is not annotated in %s genome and cloning %s seq also FAILED" \
                  " -> cannot cross-check identity with ensembl" % (gene, species, species)
        print(message)
        logging.info(message)
        return

    if not distant_seq:
        message = "\n...WARNING... : Gene %s is not annotated in %s genome " \
                  " -> cannot cross-check identity with ensembl" % (gene, species)
        print(message)
        logging.info(message)
        return

    if not seq:
        message = "\n...WARNING... : Cloning seq for gene %s from %s FAILED " \
                  " -> cannot cross-check identity with ensembl" % (gene, species)
        print(message)
        logging.info(message)
        return

    s1, s2 = aln[0].seqA, aln[0].seqB
    distant_seq_dict = {}
    for position, aa in enumerate(s1):
        distant_seq_dict[position] = aa

    seq_dict = {}
    for position, aa in enumerate(s2):
        seq_dict[position] = aa

    mismatches = 0
    for position, aa in distant_seq_dict.items():
        if seq_dict[position] != aa and seq_dict[position] != "-":
            mismatches += 1

    # report if too many mismatches occurred -> accounts for more than 1% of mismatching aa (NOT INDELS)
    score = mismatches/len(distant_seq_dict)

    if score < 0.01:
        message = "\n...NOTE... : Cloned seq for gene %s from %s species is %s percent identical to that found " \
                  "in Ensembl database" % (gene, species, 100 - (score * 100))
        print(message)
        logging.info(message)

    if score > 0.01:
        message = "\n...WARNING... : Cloned seq for gene %s from %s is only %s percent identical to that found " \
                  "in Ensembl database" % (gene, species, 100 - (score * 100))
        print(message)
        logging.info(message)


def get_most_distant_prot_seq(species, gene):
    """Finds protein sequence of the longest transcript from most distant species (e.g. rat for mouse)."""

    distant_seq = ""

    if species == "Mm":
        species = "mus musculus"
        release = 104
    if species == "Rn":
        species = "rattus norvegicus"
        release = 104
    if species == "Fc":
        species = "felis catus"
        release = 90
    if species == "Cf":
        species = "canis familiaris"
        release = 90

    logging.getLogger("pyensembl").setLevel(logging.WARNING)  # disables logging from pyensembl
    ensembl = pyensembl.EnsemblRelease(release, species)
    ensembl.download()  # this is suppose to bypass installing the release from outside python
    ensembl.index()  # this is suppose to bypass installing the release from outside python

    # check if the gene is annotated
    all_genes = ensembl.gene_names()
    if gene not in all_genes:
        return distant_seq

    # get all transcript ids
    all_transcripts_ids = ensembl.transcript_ids_of_gene_name(gene)

    all_transcripts_dict = {}
    length = 0

    # find all possible transcripts
    for t in all_transcripts_ids:
        all_transcripts_dict[t] = []

        # find gene id
        gene_id = ensembl.gene_ids_of_gene_name(gene)
        all_transcripts_dict[t].append(gene_id)

        # find contig
        contig = ensembl.locus_of_transcript_id(t).contig
        all_transcripts_dict[t].append(contig)

        # find strand
        strand = ensembl.locus_of_transcript_id(t).strand
        all_transcripts_dict[t].append(strand)

        # find transcript name
        transcript_name = ensembl.transcript_name_of_transcript_id(t)
        all_transcripts_dict[t].append(transcript_name)

        # find start of the transcript (not START codon)
        start = ensembl.locus_of_transcript_id(t).start
        all_transcripts_dict[t].append(start)

        # find end of the transcript (not STOP codon)
        end = ensembl.locus_of_transcript_id(t).end
        all_transcripts_dict[t].append(end)

        # find transcript
        biotype = "Protein coding"
        transcript = pyensembl.Transcript(t, transcript_name,
                contig, start, end, strand, biotype, gene_id, ensembl, support_level=None)

        # get protein sequence
        try:
            protein_sequence = transcript.protein_sequence
            if protein_sequence:
                length_current = len(protein_sequence)
                if length < length_current:
                    length = length_current
                    distant_seq = protein_sequence
            else:
                continue
        except ValueError:
            continue

    return distant_seq


def eliminate_frameshits_cds(wdir, ref_species, gene, raw_translated_path, raw_out_Gblocks, to_delete, aligner):
    """Finds frameshifts in cds post Gblocks at protein level and eliminates them from nucl and protein alignments."""

    raw_translated_filepath = raw_translated_path.replace(wdir, "")
    raw_out_Gblocks_filepath = raw_out_Gblocks.replace(wdir, "")
    all_seq_dict_protein = fasta_reader.alignment_file_to_dict(wdir, ref_species, raw_translated_filepath)
    all_seq_dict_cds = fasta_reader.alignment_file_to_dict(wdir, ref_species, raw_out_Gblocks_filepath)

    # make an empty dict to fill with headers and prot seq
    ref_protein = all_seq_dict_protein[">" + ref_species]

    # detect frameshifts
    for species, seq in all_seq_dict_protein.items():
        mismatches = 0

        for position, aa in enumerate(seq):
            # eliminate sequence > 10 consecutive aa are mismatched
            if mismatches > 10:
                to_delete.append(species)
                message = "\n...WARNING... : CDS for species : %s in %s protein contains a possible frameshift " \
                          "-> eliminated from alignment" % (species.replace(">", "").rstrip("\n"), gene)
                print(message)
                logging.info(message)
                break
            elif aa != ref_protein[position]:
                mismatches += 1
            # reset mismatches if not continuous
            else:
                mismatches = 0

    # make a new dict without frameshift cds
    final_dict_protein = {species: seq for (species, seq) in all_seq_dict_protein.items() if species not in to_delete}
    final_dict_cds = {species: seq for (species, seq) in all_seq_dict_cds.items() if species not in to_delete}

    # write proteins into a file
    translated_path = raw_translated_path.replace("_raw", "")
    with open(translated_path, "w") as f:
        for species, seq in final_dict_protein.items():
            f.write(species + "\n")
            f.write(seq + "\n")

    # write cds into a file
    out_Gblocks = raw_out_Gblocks.replace("aligned_" + aligner.upper() + "_Gblocks_",
                                          "PAML_template_").replace("_final_no_STOP", "")
    with open(out_Gblocks, "w") as f:
        for species, seq in final_dict_cds.items():
            f.write(species + "\n")
            f.write(seq + "\n")

    return translated_path, out_Gblocks, to_delete


def eliminate_all_insertions(gene_folder_path, out_msa):
    """Eliminates insertions in the reference sequence formed during alignment (Gblocks cannot handle them)"""

    correction = False

    alignment = AlignIO.read(gene_folder_path + "/" + out_msa, "fasta")
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]

    # check if ref species contains dashes
    if "-" in seqs[0][1]:

        positions = [p.start() for p in re.finditer("-", seqs[0][1])]

        # if yes, make a new file for a corrected alignment
        correction = True
        corrected_filename = out_msa.replace("raw", "corrected")
        file_path = gene_folder_path + "/" + corrected_filename
        with open(file_path, "w") as f:

            # loop through each sequence in the alignment
            for spec, seq in seqs:
                # generate a dict with positions as keys and bp as values (starts at 1)
                d = {}
                for i in range(len(seq)):
                    d[i+1] = seq[i]

                # eliminate these positions from given seq
                for position in positions:
                    d[position+1] = ""

                # write the corrected seqeunce into a new, corrected file
                new_seq = "".join([a for a in d.values()])
                f.write(">" + spec[2:] + "\n")
                f.write(new_seq + "\n")

            # these positions are confusingly printing index+1 (change that -> make positions start with 1)

        side_note = " \n...NOTE... : Insertions in positions in ref MSA deleted: %s" % str(positions)
        print(side_note)
        logging.info(side_note)

    # if there are no dashes in ref species sequence
    else:
        return correction, out_msa

    return correction, corrected_filename


def run_PAML(wdir, gene, PAML_path, control_file_name, codon_frequency):
    """Runs PAML by calling the control file"""

    # change working directory to given PAML gene folder
    os.chdir(PAML_path)

    # run PAML
    cml = codeml.Codeml(alignment="input.phy", tree="gene.tree", out_file="output_PAML", working_dir=PAML_path)
    cml.run(ctl_file=control_file_name, command=pyinstaller_compatibility.resource_path("codeml"))
    results = codeml.read("output_PAML")
    ns_sites = results.get("NSsites")

    M1a = ns_sites.get(1)
    M1a_lnL = M1a.get("lnL")

    M2a = ns_sites.get(2)
    M2a_lnL = M2a.get("lnL")

    M7 = ns_sites.get(7)
    M7_lnL = M7.get("lnL")
    
    M8 = ns_sites.get(8)
    M8_lnL = M8.get("lnL")
    
    # there are 2 free parameters in both comparisons
    df_M8_M7 = 2
    df_M2a_M1a = 2
    
    # determine the LRTs (Likelihood Ratio Tests) -> needs to be a positive value
    LRT1 = 2*(M2a_lnL - M1a_lnL)
    if LRT1 < 0:
        LRT1 = False
    
    LRT2 = 2*(M8_lnL - M7_lnL)
    if LRT2 < 0:
        LRT2 = False
    
    # run chi2 statistics to get p values
    if LRT1 is not False:
        M2a_M1a = cdf_chi2(df_M2a_M1a, LRT1)
    else:
        LRT1 = "(0)"
        M2a_M1a = 1
        message = "\n M1a model is more likely than M2a model."
        print(message)
        logging.info(message)
    
    if LRT2 is not False:
        M8_M7 = cdf_chi2(df_M8_M7, LRT2)
    else:
        LRT2 = "(0)"
        M8_M7 = 1
        message = "\n M7 model is more likely than M8 model."
        print(message)
        logging.info(message)
    
    message = "\n PAML LRTs (%s) for gene %s are : M2a v M1a - %s and M8 v M7 - %s" % (codon_frequency,
                                                                                     gene, str(LRT1), str(LRT2))
    print(message)
    logging.info(message)
    
    # restore working directory
    os.chdir(wdir)
    
    return M2a_M1a, M8_M7


def run_RAxML(gene, gene_folder_path, out_Gblocks):
    """Makes a gene tree using RAxML"""

    tree_name = gene + "_Tree"

    RAxML_cline = [pyinstaller_compatibility.resource_path("raxmlHPC"), "-f", "a", "-s", out_Gblocks, "-n",
                   tree_name, "-m", "GTRGAMMA", "-p", "12345", "-x", "12345", "-#", "100"]

    result = subprocess.call(RAxML_cline)
    all_tree_files = glob.glob("RAxML*")
    # move all the RAxML files to gene folder
    for file in all_tree_files:
        file_path = os.path.abspath(file)
        shutil.move(file_path, gene_folder_path)
    
    if result == 0:
        message = "\n Best gene tree was found for gene : %s " % gene
        print(message)

    # triggers "FileNotFound" exception that is handled
    else:
        message = "\n...PROBLEM... : Best gene tree was NOT found for : %s \n" % gene
        print(message)
        logging.info(message)
    
    best_tree_path = gene_folder_path + "/RAxML_bestTree." + gene + "_Tree"
    
    # remove dashes from the tree to match species names in phylip file
    with open(best_tree_path, "r") as f:
        tree = f.read()
        tree_no_dashes = tree.replace("_", "")
    
    with open(best_tree_path, "w") as f:
        f.write(tree_no_dashes)

    return best_tree_path


def get_ref_cds(wdir, gene, ref_species):
    """Takes cds for reference species from the input"""

    # open according cds fasta file
    with open(wdir + "Coding_sequences/" + gene + "_" + ref_species + "_cds.fasta", "r") as f:
        sequence = ""
        cds = f.readlines()
        for line in cds[1:]:
            sequence = sequence + line.rstrip("\n")
    return sequence


def translated_frameshift_checkpoint(wdir, seqs, gene, ref_species):
    """Cheskpoint if there is a frameshift in the translated sequence of reference species"""
    
    # take the cds of ref species post-Gblocks (always the first one)
    ref_post_Gblocks = seqs[0][1]
    # make a Seq object
    ref_post_Gblocks_seq = Seq(ref_post_Gblocks)
    # translate the Seq object
    ref_post_Gblocks_seq_object = ref_post_Gblocks_seq.translate()
    # make SeqRecord
    record_post_Gblocks = SeqRecord(ref_post_Gblocks_seq_object)

    # take the cds coming from ensembl
    ref_cds = get_ref_cds(wdir, gene, ref_species)
    # make a Seq object
    ref_cds_seq = Seq(ref_cds)
    # translate the Seq object
    ref_cds_seq_object = ref_cds_seq.translate()
    # make SeqRecord
    ref_cds_record = SeqRecord(ref_cds_seq_object)

    # perform pairwise alignment
    aln = pairwise2.align.globalxx(record_post_Gblocks.seq, ref_cds_record.seq)
    
    # make a dict storing the alignments
    d = {}
    for i in range(len(aln[0][0])):
        d[i] = ()
    
    # fill the dict with paired positions
    for position in d:
        d[position] = (aln[0][0][position], aln[0][1][position])
    
    # flag any non-synonymous substitutions that indicate frameshifts
    frameshift_positions = {}
    translated_frameshift = False
    for position, pair in d.items():
        if pair[0] != pair[1] and pair[0] != "-" and pair[1] != "-":
            # sequence index starts at 1 so need to "+1" from the python indexing
            frameshift_positions[position + 1] = pair
            translated_frameshift = True
        
    return translated_frameshift, frameshift_positions


def cloned_cds_frameshift_checkpoint(wdir, ref_species, gene, filename):
    """Compares presumptive cloned cds with reference cds and eliminates rare extreme frameshits."""

    all_seq_dict = fasta_reader.alignment_file_to_dict(wdir, ref_species, filename)

    # make an empty dict to fill with headers and cds
    ref_cds = all_seq_dict[">" + ref_species]
    to_delete = []

    # penalize mismatches (mildly) and gaps (severely)
    for species, seq in all_seq_dict.items():
        matches = 0

        for position, bp in enumerate(seq):
            if bp == "-":
                matches -= 1
                continue
            if bp == ref_cds[position]:
                matches += 1

        score = matches/len(ref_cds.replace("-", ""))
        message = "Alignment score for species : %s = %s" % (species.replace(">", ""), score)
        print(message)
        logging.info(message)

        # flag cds with rare extreme frameshifts or just missing too many exons
        # ADDED on 03/08/2023 (from 0.69 and no ref species check
        if ref_species == "Dme" and score < 0.60:
            to_delete.append(species)
            message = "...WARNING... : CDS for species : %s in %s gene is either too divergent " \
                      "or contains too many gaps -> eliminated from alignment" % (species.replace(">", "").rstrip("\n"), gene)
            print(message)
            logging.info(message)

        elif ref_species != "Dme" and score < 0.69:
            to_delete.append(species)
            message = "...WARNING... : CDS for species : %s in %s gene is either too divergent " \
                      "or contains too many gaps -> eliminated from alignment" % (species.replace(">", "").rstrip("\n"), gene)
            print(message)
            logging.info(message)

    # make a new dict without frameshift cds
    final_dict = {species: seq for (species, seq) in all_seq_dict.items() if species not in to_delete}

    # write it into a file
    with open(wdir + filename, "w") as f:
        for species, seq in final_dict.items():
            f.write(species + "\n")
            f.write(seq + "\n")

    return to_delete


def translate_Gblocks(wdir, gene_folder_path, raw_out_Gblocks_filename, gene, ref_species):
    """Translates Gblocks trimmed alignments"""
    
    filepath_to_translate = gene_folder_path + "/" + raw_out_Gblocks_filename
    # read the fasta alignment file 
    alignment = AlignIO.read(filepath_to_translate, "fasta")
    # make a list of all headers and sequences as strings
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]
    # make a filepath
    raw_translated_path = filepath_to_translate.rstrip(".fasta") + "_translated_raw.fasta"
    
    # check frameshifts in the ref species aa seq post-Gblocks
    # doesnt work if pairwise2 doesnt return alignement (Clasp1)
    translated_frameshift, frameshift_positions = translated_frameshift_checkpoint(wdir, seqs, gene, ref_species)
    if translated_frameshift is True:
        side_note = "...WARNING... : Frameshift in ref cds post-Gblocks detected: \n%s" % frameshift_positions
        print(side_note)
        logging.info(side_note)
    
    with open(raw_translated_path, "w") as f:
        # translate each sequence
        for s in seqs:            
            sequence = s[1]
            # make a Seq object
            coding_dna = Seq(sequence)
            # translate the Seq object
            seq_object = coding_dna.translate()
            # make SeqRecord
            record = SeqRecord(seq_object)
            # get a translated string using format method (only for SeqRecords)
            translated = record.format("fasta").lstrip("><unknown id> <unknown description\n")
            # write headers and sequences
            f.write(s[0].replace("_", "") + "\n")
            f.write(translated)

    # return the path to the translated alignment
    return raw_translated_path, filepath_to_translate


def convert_to_phylip(out_Gblocks):
    """Converts fasta files to phylip format"""

    phylip_path = out_Gblocks.replace(".fasta", ".phy")
    AlignIO.convert(out_Gblocks, "fasta", phylip_path, "phylip")


def run_seqret(gene, out_Gblocks):  # DEPRECATED
    """Takes a fasta squence aligment and converts it into phylip format"""

    phylip_path = out_Gblocks.replace(".fasta", ".phy")
    # PYINSTALLER: Add the pyinstaller directory to the system path so that the binary line 'which seqret' can work.
    # If not deprecated, uncomment following line to add file location to operating system path:
    # os.environ["PATH"] += os.pathsep + pyinstaller_compatibility.resource_path("")
    seqret_cline = [pyinstaller_compatibility.resource_path("seqret"),
                    "-sequence", out_Gblocks, "-osformat2", "phylipnon", "-outseq", phylip_path]
    result = subprocess.call(seqret_cline)
    if result == 0:
        message = "\n Phylip format was created for gene : %s " % gene
        print(message)
        return phylip_path
    else:
        message = "\n PROBLEM with making phylip for gene : %s " % gene
        print(message)
        logging.info(message)
      

def STOP_remover(gene_folder_path, no_dashes_out_msa, gene, aligner):
    """Removes columns (all species) following the earliest detected STOP codon"""
    
    to_trim = {}
    max_trim_position = 30
    c_term_bp_to_trim = 0
    path_to_trim = gene_folder_path + "/" + no_dashes_out_msa
    with open(path_to_trim, "r") as f:
        file = f.readlines()
        
        for line in file:
            
            if line.startswith(">"):
                seq = ""
                header = line.rstrip("\n")
                to_trim[header] = seq
            
            else:
                s = line.rstrip("\n")
                seq = seq + s
                to_trim[header] = seq
    
    for species, cds in to_trim.items():
        
        # split the end of each cds to look for early STOP codon
        c_term = cds[-max_trim_position:]
        c_term_codons = [c_term[i:i+3] for i in range(0, len(c_term), 3)]
            
        for codon in c_term_codons:
            if codon.upper() in ["TAA", "TGA", "TAG"]:
                # find how many bp to remove from c_term
                new_c_term_bp_to_trim = (len(c_term_codons) - c_term_codons.index(codon)) * 3
                
                # make sure that the earliest STOP is the final trim position
                if new_c_term_bp_to_trim > c_term_bp_to_trim:
                    c_term_bp_to_trim = new_c_term_bp_to_trim
    
    final_cds_no_STOP = {}
    for species, cds in to_trim.items():
        # subtract number of dashes to remove longer more bases (if STOP detected)
        if c_term_bp_to_trim != 0:
            final_cds_no_STOP[species] = cds[:-c_term_bp_to_trim]
        # no STOP codon detected in any species (ex. C-term tip missing; CENPX last exon (5) is microexon in primates)
        else:
            final_cds_no_STOP[species] = cds

    final_cds_file_no_STOP = "aligned_" + aligner.upper() + "_" + gene + "_final_no_STOP.fasta"
    with open(final_cds_file_no_STOP, "w") as f:
    
        for species, cds in final_cds_no_STOP.items():
            f.write(species + "\n")
            f.write(cds + "\n")
    
    return final_cds_file_no_STOP


def post_Gblocks_STOP_remover(gene_folder_path, raw_out_Gblocks_filename, gene):
    """Removes a single column from all species in case a STOP is detected post-Gblocks (rare).
    Forces conserved alignment and prevents frameshifts"""
    
    post_Gblocks_path = gene_folder_path + "/" + raw_out_Gblocks_filename
    
    all_cds = {}
    all_cds_no_STOP = {}
    premature_STOP_codon_positions = set()  # doesnt allow repetitions
    species_with_STOP = []
    
    # read species and cds into a dict to ease search
    with open(post_Gblocks_path, "r") as f:
        file = f.readlines()
        
        for line in file:
            
            if line.startswith(">"):
                seq = ""
                header = line.rstrip("\n")
                all_cds[header] = seq
            
            else:
                s = (line.rstrip("\n").replace(" ", "")).upper()
                seq = seq + s
                all_cds[header] = seq
    
    # chop cds into codons (all are in frame) and mark STOPs
    for header, cds in all_cds.items():

        codons_list = [cds[i:i+3] for i in range(0, len(cds), 3)]
        STOP_positions = [position for position, codon in enumerate(codons_list) if codon in ["TAA", "TGA", "TAG"]]

        if STOP_positions:
            
            species = header.replace(">", "")
            species_with_STOP.append(species)
            
            for position in STOP_positions:
                premature_STOP_codon_positions.add(position)
    
    # edit all coding sequences removing STOP columns
    for header, cds in all_cds.items():
        
        all_codons = [cds[i:i+3] for i in range(0, len(cds), 3)]
        cds_no_STOP = ""
        
        # join all the codon except STOP positions
        all_cds_no_STOP[header] = "".join([cds_no_STOP + codon for position, codon in enumerate(all_codons) \
                                   if position not in premature_STOP_codon_positions])
    
    # overwrite the previous out_Gblocks file
    with open(post_Gblocks_path, "w") as f:
        
        for header, cds in all_cds_no_STOP.items():
            f.write(header + "\n")
            f.write(cds + "\n")
    
    # issue a warning if STOPs were present
    if species_with_STOP:

        STOP_positions = sorted([(position * 3) + 1 for position in premature_STOP_codon_positions])
        
        message = "\n !!! WARNING !!! - CDS for: %s contains early STOP codons " \
                  "starting at nucleotide position in Gblocks_final_no_STOP alignment: %s in species: %s \n         " \
                  "----> they were removed in each species forcing conserved alignment" \
                  % (gene, STOP_positions, species_with_STOP)
        print(message)
        logging.info(message)
    
    return post_Gblocks_path


def run_Gblocks(final_cds_file_no_STOP, gene, result_path, aligner):
    """Runs Gblocks software which removes indels in frame from the final alignment"""

    in_filepath = result_path + gene + "/" + final_cds_file_no_STOP
    raw_out_Gblocks_filename = "aligned_" + aligner.upper() + "_Gblocks_" + gene + "_final_no_STOP.fasta"
    # run Gblocks with options: codon ("-t=c") and dont save html file ("-p=n")
    Gblocks_cline = [pyinstaller_compatibility.resource_path("Gblocks"), in_filepath, "-t=c", "-p=n"]
    result = subprocess.run(Gblocks_cline, capture_output=True)
    # get stdout -> decode to string from bit -> split by lines and get the 5th and 6th
    message = (result.stdout.decode('utf-8').split("\n"))[5:7]
    # log stdout
    print(message)
    logging.info(message)
    # add fasta extension
    try:
        os.rename(in_filepath + "-gb", raw_out_Gblocks_filename)
    except FileNotFoundError:
        message = "\n...WARNING... : Failed PAML analysis for : %s -> probably not enough " \
                  "species in the alignment (e.g. poor alignment of repetitive regions)" % gene
        print(message)
        logging.info(message)
        return
    
    # returns the filename after Gblocks
    return raw_out_Gblocks_filename

              
def align_final_cds(wdir, gene, result_path):
    """Aligns cloned coding sequences that pass coverage threshold (default 90 percent)"""

    gene_folder_path = result_path + gene + "/"
    in_filepath = gene_folder_path + gene + "_final.fasta"
    out_msa = gene + "_raw_nucleotide_alignment.fasta"

    aligner = "mafft"

    # define which aligner is used
    if aligner == "mafft":   # runs mafft in G-INS-i with VSM (Katoh et al., 2016 Bioinformatics) added 07/10/2022
        cmd = [pyinstaller_compatibility.resource_path("mafft-ginsi"), "--allowshift", "--unalignlevel", "0.8",
               "--out", gene_folder_path + out_msa, in_filepath]
        subprocess.call(cmd)
        # returns the filename after MSA
        return out_msa, aligner

    if aligner == "guidance":  # NOT SUPPORTED YET
        guidance_path = "/Users/damian/PycharmProjects/Freeda_pyinstaller_04_29_2022/include_mac/www/Guidance/guidance.pl"
        cmd = ["perl", guidance_path, "--seqFile", in_filepath, "--msaProgram", "MAFFT", "--mafft",
               pyinstaller_compatibility.resource_path("mafft-ginsi"), "--seqType", "nuc",
               "--outDir", gene_folder_path + out_msa, "--outOrder", "as_input"]
        subprocess.call(cmd)
        # returns the filename after MSA
        return out_msa, aligner

    elif aligner == "muscle":  # POOR ALIGNMENT USUALLY -> DEPRECATED

        cmd = ['muscle', "-in", in_filepath, "-quiet", "-maxiters", "2", "-out", gene_folder_path + out_msa]
        subprocess.call(cmd)
        # need to reorder seqs post msa
        fasta_reader.reorder_alignment(in_filepath, gene_folder_path + out_msa)

        # returns the filename after MSA
        return out_msa, aligner

    elif aligner == "prank":  # TAKES TOO LONG

        cline = PrankCommandline(d=in_filepath,
                                       o=out_msa,  # prefix only!
                                       f=8)  # FASTA output
        # record standard output and standard error
        stdout, stderr = cline()
        shutil.move(wdir + out_msa + ".best.fas", gene_folder_path + out_msa)
        fasta_reader.reorder_alignment(in_filepath, gene_folder_path + out_msa)

        # returns the filename after MSA
        return out_msa.replace(".best.fas", ""), aligner


