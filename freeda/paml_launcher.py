#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 10:51:27 2021

@author: damian
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 21:48:45 2021

@author: damian

Analyses the final cds, gets a gene tree based on translated cds and runs PAML.

"""

from freeda import fasta_reader
from freeda import genomes_preprocessing
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import datetime
import time
import glob
import os
import re
import logging
import shutil
import subprocess
import copy


def analyse_final_cds(wdir, ref_species, result_path, all_proteins):

    failed_paml = []
    start_time = time.time()
    day = datetime.datetime.now().strftime("-%m-%d-%Y-%H-%M")
        
    # initiate log file to record PAML analysis by reseting the handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    PAML_logfile_name = "PAML" + day + ".log"
    logging.basicConfig(filename=PAML_logfile_name, level=logging.INFO, format="%(message)s")
    
    # make empty dataframe to store PAML results
    #PAML_df = pd.DataFrame({"Protein Name" : \
    #                        "p-value" : \
    #                        "M2a vs M1a"})
    
    # make an empty template for all possible species
    all_species = [names[0] for names in genomes_preprocessing.get_names(ref_species)]

    all_species_dict = {ref_species: ""}
    for species in all_species:
        all_species_dict[species] = ""

    #with open("species.txt", "r") as f:
    #    all_species[ref_species] = ""
    #    for species in f.readlines():
    #        all_species[species.rstrip("\n")] = ""
    
    nr_of_species_total_dict = {}
    prots_under_pos_sel = []

    for protein in all_proteins:
        
        # check if this protein was already analysed
        if os.path.isdir(result_path + protein + "/" + "PAML_" + protein):
           message = "\n################\n\n PAML analysis has been already performed for : %s (skipping)" % protein
           print(message)
           logging.info(message)
           
           # get how many species were analysed
           path_to_final_species = result_path + protein + "/" + protein + "_final.fasta"
           with open(path_to_final_species, "r") as f:
               nr_of_species = f.read().count(">")
               nr_of_species_total_dict[protein] = nr_of_species
           
        # otherwise proceed with the analysis
        else:
            # need to use deepcopy function to make an actual dictionary copy
            final_species = copy.deepcopy(all_species_dict)
            final_species_headers = []
        
            # read each cds fasta file and put them into the empty final_species dict
            cds_file_path = result_path + protein + ".fasta"
        
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
                        ref_head = line.split(protein + "_")[1].rstrip("\n")
                        continue
                
                    # second line is always the ref species sequence
                    if line_nr == 2:
                        # remove STOP codon if present
                        # ref_seq = STOP_remover(line.rstrip("\n"))
                        ref_seq = line.rstrip("\n")
                        final_species[ref_head] = ref_seq
                        continue
                
                    # next lines are headers and sequences of cloned cds
                    if line.startswith(">"):
                        head = line.split(protein + "_")[1].rstrip("\n")
                        continue
                
                    else:
                        # remove STOP codon if present
                        # seq = STOP_remover(line.rstrip("\n"))
                        seq = line.rstrip("\n")
                        seq_no_dashes = seq.replace("-", "")
                        if len(seq_no_dashes) / len(ref_seq) >= 0.90:
                            final_species[head] = seq
                            final_species_headers.append(head)
               
            # log species
            final_species_headers.insert(0, ref_species)
            nr_of_species = len(final_species_headers)
            nr_of_species_total_dict[protein] = nr_of_species

            message = "\n\n --------- * %s * --------- \n\n" % protein
            print(message)
            logging.info(message)
        
            # write all sequences into a file by reading through the dictionary
            final_cds_file = protein + "_final.fasta"
            with open(final_cds_file, "w") as f:
                for species in final_species:
                    
                    if final_species[species] == "":
                        continue
                    else:
                        f.write(">" + species + "\n")
                        f.write(final_species[species] + "\n")
        
        
            protein_folder_path = result_path + protein
            
            shutil.move(cds_file_path, protein_folder_path)
            shutil.move(final_cds_file, protein_folder_path)
        
            # generate a PAML folder for a given protein
            PAML_path = protein_folder_path + "/PAML_" + protein
            os.makedirs(PAML_path)
        
            # copy control_file from working directory into protein_folder_path
            control_file = "control_file.ctl"
            shutil.copy("control_file.ctl", PAML_path + "/" + control_file)
        
            # align the final cds sequences
            out_mafft = align_final_cds(protein, final_cds_file, result_path)

            # check for rare frameshifts in cloned cds
            cloned_cds_frameshift_checkpoint(wdir, ref_species, protein, out_mafft)
            shutil.move(out_mafft, protein_folder_path)

            # check and eliminate insertions that cause dashes in ref species
            correction, corrected_filename = eliminate_all_insertions(protein_folder_path, out_mafft)
            if correction is True:
                side_note = " \n---- Insertions detected in " + protein + " alignment -> " \
                    "these bp positions were removed in all species forcing conserved alignment\n"
                print(side_note)
                logging.info(side_note)
                no_dashes_out_mafft = corrected_filename
            else:
                no_dashes_out_mafft = out_mafft
        
            # remove STOP codons
            final_cds_file_no_STOP = STOP_remover(protein_folder_path, no_dashes_out_mafft, protein)
            shutil.move(final_cds_file_no_STOP, protein_folder_path)
        
            # run gBLOCK
            raw_out_Gblocks_filename = run_Gblocks(final_cds_file_no_STOP, protein, result_path)
            # Gblocks will fail if not enough species in alignment -> go to next protein
            if raw_out_Gblocks_filename is None:
                #nr_of_species_total_dict = False
                #PAML_logfile_name = False
                #day = False
                failed_paml.append(protein)
                #prots_under_pos_sel = False
                #return nr_of_species_total_dict, PAML_logfile_name, day, failed_paml, prots_under_pos_sel
                continue

            shutil.move(wdir + raw_out_Gblocks_filename, protein_folder_path)

            # double check for artificial STOP codons introduced by Gblocks -> force conserved alignment
            # it overwrites the previous out_Gblocks file (same name)
            post_Gblocks_STOP_remover(protein_folder_path, raw_out_Gblocks_filename, protein)

            # translate all sequences (inside of the protein folder)
            raw_translated_path, filepath_to_translate = translate_Gblocks(wdir, protein_folder_path,
                                                                        raw_out_Gblocks_filename, protein, ref_species)

            # check for frameshifts in protein seq and eliminate species from protein and cds
            translated_path, out_Gblocks, to_delete = eliminate_frameshits_cds(wdir, ref_species, protein,
                                                                    raw_translated_path, filepath_to_translate)

            # correct number of species in the final analysis
            nr_of_species = nr_of_species - len(to_delete)
            final_species_headers = [species for species in final_species_headers if species not in to_delete]
            message = "\n Final species cloned and aligned (+ ref) for %s : %s %s \n" \
                      % (protein, str(nr_of_species), str(final_species_headers))
            print(message)
            logging.info(message)

            # run seqret (file is already in protein folder)
            phylip_path = run_seqret(protein, out_Gblocks)
            shutil.copy(phylip_path, PAML_path + "/input.phy")
            
            # run RAxML (and move all the RAxML files to protein folder)
            try:
                best_tree_path = run_RAxML(protein, protein_folder_path, translated_path)

            except FileNotFoundError:
                failed_paml.append(protein)
                message = "...PROBLEM... : The best gene tree was NOT found for %s " \
                          "-> likely too few cds pass coverage threshold (90 percent default)" % protein
                print(message)
                logging.info(message)
                continue

            shutil.copy(best_tree_path, PAML_path + "/gene.tree")

            # rename and copy the final protein alignment into results
            shutil.copy(translated_path, result_path + protein + "_protein_alignment.fasta")
        
            # run PAML
            message = "\n.........Running PAML for protein: %s.........\n" % protein
            print(message)
            
            M2a_M1a, M8_M7 = run_PAML(wdir, protein, PAML_path, control_file)

            if M8_M7 is not None and M8_M7 < 0.05:
                prots_under_pos_sel.append(protein)

            message = "\n -> PAML p-values for protein %s : M2a v M1a - %s and M8 v M7 - %s" \
                % (protein, str(M2a_M1a), str(M8_M7))
            print(message)
            logging.info(message)

    shutil.move(wdir + PAML_logfile_name, result_path)

    # for now final analysis complete message logs into the PAML log file -> change that later
    
    # mark the end of the analysis
    message = ("\n --------------->  PAML analysis completed in %s minutes or %s hours"
               % ((time.time() - start_time)/60, (time.time() - start_time)/60/60))
    print(message)
    logging.info(message)
    
    return nr_of_species_total_dict, PAML_logfile_name, day, failed_paml, prots_under_pos_sel


def eliminate_frameshits_cds(wdir, ref_species, protein, raw_translated_path, raw_out_Gblocks):
    """Finds frameshifts in cds post Gblocks at protein level and eliminates them from nucl and protein alignments."""

    raw_translated_filepath = raw_translated_path.replace(wdir, "")
    raw_out_Gblocks_filepath = raw_out_Gblocks.replace(wdir, "")
    all_seq_dict_protein = fasta_reader.alignment_file_to_dict(wdir, ref_species, raw_translated_filepath)
    all_seq_dict_cds = fasta_reader.alignment_file_to_dict(wdir, ref_species, raw_out_Gblocks_filepath)

    # make an empty dict to fill with headers and prot seq
    ref_protein = all_seq_dict_protein[">" + ref_species]
    to_delete = []

    # detect frameshifts
    for species, seq in all_seq_dict_protein.items():
        mismatches = 0

        for position, aa in enumerate(seq):
            # eliminate sequence > 10 consecutive aa are mismatched (threshold at 10 legit mismatches in Izumo1 Rn)
            if mismatches > 10:
                to_delete.append(species)
                message = "\n...WARNING... : CDS for species : %s in %s protein contains a frameshift " \
                          "-> eliminated from alignment" % (species.rstrip("\n"), protein)
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
    out_Gblocks = raw_out_Gblocks.replace("aligned_MAFFT_Gblocks_", "PAML_template_").replace("_final_no_STOP", "")
    with open(out_Gblocks, "w") as f:
        for species, seq in final_dict_cds.items():
            f.write(species + "\n")
            f.write(seq + "\n")


    return translated_path, out_Gblocks, to_delete


def eliminate_all_insertions(protein_folder_path, out_mafft):

    correction = False
    
    alignment = AlignIO.read(protein_folder_path + "/" + out_mafft, "fasta")
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]
    
    # check if ref species contains dashes
    if "-" in seqs[0][1]:
        
        positions = [p.start() for p in re.finditer("-", seqs[0][1])]
        
        # if yes, make a new file for a corrected alignment
        correction = True
        corrected_filename = "corrected_" + out_mafft
        file_path = protein_folder_path + "/" + corrected_filename
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
            
        side_note = " \n WARNING : Insertions in positions in ref MSA deleted: %s" % str(positions)
        print(side_note)
        logging.info(side_note)
                    
    # if there are no dashes in ref species sequence
    else:
        return correction, out_mafft
    
    return correction, corrected_filename


# MODIFY THE INSERTION FUNCTION TO ELIMINATE ALSO 3% != 0 insertions > 1 (mostly artificial insertions)



def run_PAML(wdir, protein, PAML_path, control_file):

    from Bio.Phylo.PAML import codeml
    from Bio.Phylo.PAML.chi2 import cdf_chi2
    
    # change working directory to given PAML protein folder
    os.chdir(PAML_path)
    
    # run PAML
    cml = codeml.Codeml(alignment="input.phy", tree="gene.tree", out_file="output_PAML", working_dir=PAML_path)
    cml.run(control_file)
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
    if LRT1 != False:
        M2a_M1a = cdf_chi2(df_M2a_M1a, LRT1)
    else:
        LRT1 = None
        M2a_M1a = None
        message = "\n M1a model is more likely than M2a model."
        print(message)
        logging.info(message)
    
    if LRT2 != False:
        M8_M7 = cdf_chi2(df_M8_M7, LRT2)
    else:
        LRT2 = None
        M8_M7 = None
        message = "\n M7 model is more likely than M8 model."
        print(message)
        logging.info(message)
    
    message = "\n PAML LRTs for protein %s are : M2a v M1a - %s and M8 v M7 - %s" \
            % (protein, str(LRT1), str(LRT2))
    print(message)
    logging.info(message)
    
    # restore working directory
    os.chdir(wdir)
    
    return M2a_M1a, M8_M7
    

def run_RAxML(protein, protein_folder_path, translated_path):
    tree_name = protein + "_Tree"
    RAxML_cline = ["raxmlHPC", "-f", "a", "-s", translated_path, "-n", tree_name,
                   "-m", "PROTGAMMAAUTO", "-p", "1985", "-x", "2020", "-#", "100"]
    result = subprocess.call(RAxML_cline)
    all_tree_files = glob.glob("RAxML*")
    # move all the RAxML files to protein folder
    for file in all_tree_files:
        file_path = os.path.abspath(file)
        shutil.move(file_path, protein_folder_path)
    
    if result == 0:
        message = "\n Best gene tree was found for protein : %s " % protein
        print(message)
        #logging.info(message)

    # triggers "FileNotFound" exception that is handled
    else:
        message = "\n...PROBLEM... : Best gene tree was NOT found for : %s \n" % protein
        print(message)
        logging.info(message)
    
    best_tree_path = protein_folder_path + "/RAxML_bestTree." + protein + "_Tree"
    
    # remove dashes from the tree to match species names in phylip file
    with open(best_tree_path, "r") as f:
        tree = f.read()
        tree_no_dashes = tree.replace("_", "")
    
    with open(best_tree_path, "w") as f:
        f.write(tree_no_dashes)

    return best_tree_path


def get_ref_cds(wdir, protein_name, ref_species):

    # open according cds fasta file
    with open(wdir + "Coding_sequences/" + protein_name + "_" + ref_species + "_cds.fasta", "r") as f:
        sequence = ""
        cds = f.readlines()
        for line in cds[1:]:
            sequence = sequence + line.rstrip("\n")
    return sequence


def translated_frameshift_checkpoint(wdir, seqs, protein, ref_species):
    from Bio import pairwise2
    
    # take the cds of ref species post-Gblocks (always the first one)
    ref_post_Gblocks = seqs[0][1]
    # make a Seq object
    ref_post_Gblocks_seq = Seq(ref_post_Gblocks)
    # translate the Seq object
    ref_post_Gblocks_seq_object = ref_post_Gblocks_seq.translate()
    # make SeqRecord
    record_post_Gblocks = SeqRecord(ref_post_Gblocks_seq_object)

    # take the cds coming from ensembl
    ref_cds = get_ref_cds(wdir, protein, ref_species)
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


def cloned_cds_frameshift_checkpoint(wdir, ref_species, protein, out_mafft):
    """Compares presumptive cloned cds with reference cds and eliminates rare extreme frameshits."""

    all_seq_dict = fasta_reader.alignment_file_to_dict(wdir, ref_species, out_mafft)

    # make an empty dict to fill with headers and cds
    ref_cds = all_seq_dict[">" + ref_species]
    to_delete = []

    # penalize mismatches and gaps
    for species, seq in all_seq_dict.items():
        matches = 0

        for position, bp in enumerate(seq):
            if bp == "-":
                matches -= 1
                continue
            if bp == ref_cds[position]:
                matches += 1

        score = matches/len(ref_cds.replace("-", ""))
        message = "Alignment score for species : %s = %s" % (species.rstrip("\n"), score)
        print(message)
        logging.info(message)

        # flag cds with rare extreme frameshifts
        if score < 0.55: # threshold set by legi >Gd = 0.5879396984924623 Izumo1
            to_delete.append(species)
            message = "\n...WARNING... : CDS for species : %s in %s protein contains a frameshift " \
                      "-> eliminated from alignment" % (species.rstrip("\n"), protein)
            print(message)
            logging.info(message)

    # make a new dict without frameshift cds
    final_dict = {species: seq for (species, seq) in all_seq_dict.items() if species not in to_delete}

    # write it into a file
    with open(wdir + out_mafft, "w") as f:
        for species, seq in final_dict.items():
            f.write(species + "\n")
            f.write(seq + "\n")


def translate_Gblocks(wdir, protein_folder_path, raw_out_Gblocks_filename, protein, ref_species):
    
    filepath_to_translate = protein_folder_path + "/" + raw_out_Gblocks_filename
    # read the fasta alignment file 
    alignment = AlignIO.read(filepath_to_translate, "fasta")
    # make a list of all headers and sequences as strings
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]
    # make a filepath
    raw_translated_path = filepath_to_translate.rstrip(".fasta") + "_translated_raw.fasta"
    
    # check frameshifts in the ref species aa seq post-Gblocks
    # doesnt work if pairwise2 doesnt return alignement (Clasp1)
    translated_frameshift, frameshift_positions = translated_frameshift_checkpoint(wdir, seqs, protein, ref_species)
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


def run_seqret(protein, out_Gblocks):
    #in_filepath = protein_folder_path + "/" + out_Gblocks
    phylip_path = out_Gblocks.replace(".fasta", ".phy")
    seqret_cline = ["seqret", "-sequence", out_Gblocks,
                    "-osformat2", "phylipnon", "-outseq", phylip_path]
    result = subprocess.call(seqret_cline)
    if result == 0:
        message = "\n Phylip format was created for protein : %s " % protein
        print(message)
        #logging.info(message)
        return phylip_path
    else:
        message = "\n PROBLEM with making phylip for protein : %s " % protein
        print(message)
        logging.info(message)
      

def STOP_remover(protein_folder_path, no_dashes_out_mafft, protein_name):
    
    to_trim = {}
    max_trim_position = 30
    c_term_bp_to_trim = 0
    path_to_trim = protein_folder_path + "/" + no_dashes_out_mafft
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

    final_cds_file_no_STOP = "aligned_MAFFT_" + protein_name + "_final_no_STOP.fasta"
    with open(final_cds_file_no_STOP, "w") as f:
    
        for species, cds in final_cds_no_STOP.items():
            f.write(species + "\n")
            f.write(cds + "\n")
    
    return final_cds_file_no_STOP


def post_Gblocks_STOP_remover(protein_folder_path, raw_out_Gblocks_filename, protein):
    # checks if artificial STOP was formed after Gblocks
    
    post_Gblocks_path = protein_folder_path + "/" + raw_out_Gblocks_filename
    
    all_cds = {}
    all_cds_no_STOP = {}
    artificial_STOP_codon_positions = set() # doesnt allow repetitions
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
                artificial_STOP_codon_positions.add(position)
    
    # edit all coding sequences removing STOP columns
    for header, cds in all_cds.items():
        
        all_codons = [cds[i:i+3] for i in range(0, len(cds), 3)]
        cds_no_STOP = ""
        
        # join all the codon except STOP positions
        all_cds_no_STOP[header] = "".join([cds_no_STOP + codon for position, codon in enumerate(all_codons) \
                                   if position not in artificial_STOP_codon_positions])
    
    # overwrite the previous out_Gblocks file
    with open(post_Gblocks_path, "w") as f:
        
        for header, cds in all_cds_no_STOP.items():
            f.write(header + "\n")
            f.write(cds + "\n")
    
    # issue a warning if STOPs were present
    if species_with_STOP:
        
        STOP_positions = sorted([(position * 3) + 1 for position in artificial_STOP_codon_positions])
        
        message = "\n !!! WARNING !!! - CDS for: %s contains early STOP codons " \
        "starting at nucleotide position in Gblocks_final_no_STOP alignment: %s in species: %s \n           " \
        "----> they were removed in each species forcing conserved alignment" \
                % (protein, STOP_positions, species_with_STOP) 
        print(message)
        logging.info(message)
    
    return post_Gblocks_path


def run_Gblocks(final_cds_file_no_STOP, protein, result_path):

    in_filepath = result_path + protein + "/" + final_cds_file_no_STOP
    raw_out_Gblocks_filename = "aligned_MAFFT_Gblocks_" + protein + "_final_no_STOP.fasta"
    # run Gblocks with options: codon ("-t=c") and dont save html file ("-p=n")
    Gblocks_cline = ["Gblocks", in_filepath, "-t=c", "-p=n"]
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
                  "species in the alignment (e.g. poor alignment of repetitive regions)" % protein
        print(message)
        logging.info(message)
        return
    
    # returns the filename after Gblocks
    return raw_out_Gblocks_filename

              
def align_final_cds(protein, final_cds_file, result_path):
    
    in_filepath = result_path + protein + "/" + protein + "_final.fasta"
    out_mafft = "aligned_MAFFT_" + final_cds_file
    # run mafft
    mafft_cline = MafftCommandline(input=in_filepath)
    # record standard output and standard error
    stdout, stderr = mafft_cline()
    # make a post-MSA file using out_filename
    with open(out_mafft, "w") as f:
        f.write(stdout)
    
    # returns the filename after MAFFT
    return out_mafft


"""


    all_seq_dict = {}
    cds_seq = ""

    # write mafft output into dict
    with open(wdir + out_mafft, "r") as f:
        file = f.readlines()

        for line in file:

            if line.startswith(">" + ref_species):
                header = line.rstrip("\n")
                continue

            if not line.startswith(">"):
                cds_seq += line.rstrip("\n")
                continue

            # adds ref species header the first time its executed
            if line.startswith(">"):
                all_seq_dict[header] = cds_seq
                header = line.rstrip("\n")
                cds_seq = ""
                continue

        # record the final species
        all_seq_dict[header] = cds_seq


def remove_non_ACGT_codons(result_path, out_Gblocks):
    
    non_ACGT = ["N","Y","R","W","S","K","M","D","H","V","B","X"]
    cleaned_Gblocks = open(out_Gblocks[:-6] + "_cleaned.fasta") 
    out_Gblocks = result_path + "aligned_MAFFT_Gblocks_Terf1_final.fasta"
    with open(out_Gblocks, "r") as f:
        line_nr = 0
        file = f.readlines()
        for line in file:
            line_nr += 1
            if line.startswith('>'):
                pass
                #cleaned_Gblocks.write(line)
            elif any(n in line.upper() for n in non_ACGT):
                print(line)
                print(line_nr)
            else:
                pass
                #cleaned_Gblocks.write(line)
"""



