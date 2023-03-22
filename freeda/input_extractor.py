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

Extracts input using pyensembl API that communicates with Ensembl database

"""

from freeda import tblastn, genomes_preprocessing, pyinstaller_compatibility, fasta_reader
from Bio import SeqIO
import pyensembl
import pybedtools
import os
import shutil
import logging
import requests
from requests.exceptions import HTTPError


# PYINSTALLER: Set bedtools path to a bedtools folder in the FREEDA directory.
if pyinstaller_compatibility.is_bundled():
    pybedtools.helpers.set_bedtools_path(pyinstaller_compatibility.resource_path("bedtools/bin"))
else:
    pybedtools.helpers.set_bedtools_path("")

rules = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
         "Y": "R", "R": "Y", "W": "W", "S": "S", "K": "M", "M": "K",
         "D": "H", "H": "D", "V": "B", "B": "V", "X": "X"}


def correct_for_microexons(wdir, ref_species, gene, microexons, missing_bp_list, transcript):
    """Corrects the automatic input - exons and cds. IT CANNOT HANDLE TWO MICROEXONS ONE AFTER THE OTHER"""

    input_correct = True

    path_to_exons = wdir + "Exons/"
    path_to_cds = wdir + "Coding_sequences/"

    # these exons dont have microexons anymore (expected_exons not used here)
    ref_exons, expected_exons = fasta_reader.get_ref_exons(wdir, gene, ref_species, at_input=True)

    with open(path_to_exons + gene + "_" + ref_species + "_exons.fasta", "w") as f:

        # make a shallow copy (doesnt modify the original)
        final_exons = ref_exons.copy()

        # update final_exons dict for each microexon
        for microexon in microexons:

            bp_count = 0
            microexon_nr = microexon[0]
            microexon_length = missing_bp_list.pop(0)

            # keeps track of the frame before and after microexon
            trimmed_bp = 0

            if microexon_nr == 1:
                bp_count += microexon_length

            for exon in final_exons:

                # final_exons dict does not have microexons !!!
                ex_header, ex_seq, ex_length = final_exons[exon]

                # find the exon preceding a microexon
                if exon == microexon_nr - 1:

                    # need to add it to bp_count first before trimming
                    bp_count += ex_length

                    # first frame -> no need for corrections
                    if bp_count % 3 == 0:
                        pass

                    # second frame -> need to trim one bp
                    elif bp_count % 3 == 1:
                        bp_count -= 1
                        trimmed_bp = 1
                        # update exon length
                        ex_seq = ex_seq[:-1]
                        ex_length = ex_length - 1

                    # third frame -> need to trim two bps
                    else:
                        bp_count -= 2
                        trimmed_bp = 2
                        # update exon length
                        ex_seq = ex_seq[:-2]
                        ex_length = ex_length - 2

                    # update ref exons dictionary (its a shallow copy)
                    final_exons[exon] = ex_header, ex_seq, ex_length

                    continue

                # find the exon proceeding a microexon
                if exon == microexon_nr + 1:

                    # removing microexon does not alter the frame
                    if (microexon_length + trimmed_bp) % 3 == 0:
                        bp_count += ex_length

                    # removing microexon alters the frame
                    elif (microexon_length + trimmed_bp) % 3 == 1:
                        ex_seq = ex_seq[2:]
                        bp_count += ex_length - 2
                        # update exon length
                        ex_length = ex_length - 2

                    # removing microexon alters the frame
                    else:
                        ex_seq = ex_seq[1:]
                        bp_count += ex_length - 1
                        # update exon length
                        ex_length = ex_length - 1

                    # update ref exons dictionary (its a shallow copy)
                    final_exons[exon] = ex_header, ex_seq, ex_length

                    continue

                else:
                    bp_count += ex_length

                    # update ref exons dictionary (its a shallow copy)
                    final_exons[exon] = ex_header, ex_seq, ex_length

        corrected_cds = ""
        for exon in final_exons:
            f.write(final_exons[exon][0] + "\n")
            f.write(final_exons[exon][1] + "\n")
            corrected_cds += final_exons[exon][1]

    # correct cds (but not the gene -> better to get micrexon matches to extend contig)
    with open(path_to_cds + gene + "_" + ref_species + "_cds.fasta", "w") as f:
        f.write(">" + transcript.name + "_" + ref_species + "_cds\n")
        f.write(corrected_cds + "\n")

    if len(corrected_cds) % 3 != 0:
        input_correct = False
        message = "\n...WARNING... : CDS of %s in ref species is NOT in frame." % gene
        logging.info(message)

        return input_correct

    return input_correct


def get_uniprot_id(ref_species, gene, ensembl):
    """Retrieves all possible uniprot ids to be matched against structure prediction from AlphaFold"""

    if ref_species == "Mm":
        ref_species_number = "10090"
    elif ref_species == "Rn":
        ref_species_number = "10116"
    elif ref_species == "Hs":
        ref_species_number = "9606"
    elif ref_species == "Fc":
        ref_species_number = "9685"
    elif ref_species == "Cf":
        ref_species_number = "9615"
    elif ref_species == "Gg":
        ref_species_number = "9031"
    elif ref_species == "Dme":
        ref_species_number = "7227"
        # need to convert the FlyBase ID to gene name for Uniprot requests
        gene = ensembl.gene_by_id(gene).gene_name

    possible_uniprot_ids = []
    # AlphaFold flags "reviewed" protein with "GENENAME_MOUSE" so name count helps to find canonical version
    name_count = 0
    best_id = None

    query = "https://rest.uniprot.org/uniprotkb/search?query=gene_exact:" + gene + "+AND+taxonomy_id:" + \
            ref_species_number + "&format=tsv"

    # get possible uniprot ids
    try:  # need to use exception handling as sometimes error in connecting to url
        data = requests.get(query)
        data.raise_for_status()  # Raises an HTTPError if the return code is 4xx or 5xx
        data_text = data.text  # gets a string from tab format
        data_list = data_text.split("\n")

        # record reviewed entry
        for element in data_list:
            # it might help to resolve multiple proteins under same alternative gene name
            current_name_count = element.count(gene) + element.count(gene.upper())
            # uniprot ID is always under "Entry" which is the first element
            possible_id = element.split("\t")[0]
            # define if entry is reviewed
            try:
                reviewed = element.split("\t")[2]
            except IndexError:  # when empty element
                continue

            if current_name_count > name_count:
                name_count = current_name_count
                best_id = possible_id
            if possible_id != "Entry" and reviewed.upper() == "REVIEWED":  # only consider "reviewed" entries
                possible_uniprot_ids.append(possible_id)

        if not possible_uniprot_ids:
            # get the first non-reviewed entry (limited chance that this entry will have an AlphaFold prediction)
            for element in data_list:
                # uniprot ID is always under "Entry" which is the first element
                possible_id = element.split("\t")[0]
                if possible_id != "Entry" and possible_id != "":
                    possible_uniprot_ids.append(possible_id)
                    break
            return possible_uniprot_ids

        # there was no entries, including non-reviewed ones
        if not possible_uniprot_ids:
            return possible_uniprot_ids

        # put the id with highest probability (name_count) last -> make a dict (remove duplicates) from list and back
        # to list, then append best_id
        id_list = list(dict.fromkeys(possible_uniprot_ids))
        # need to get rid of the last id to prevent duplications
        # added "in possible_uniprot_ids" to avoid ValueError if best_id is not in possible_uniprot_ids
        # cose unreviewed entry had more name counts for some reason (e.g. Magi1 in rodents)
        if best_id in possible_uniprot_ids:
            id_list.remove(best_id)
            id_list.append(best_id)
        possible_uniprot_ids = id_list

    except HTTPError:
        pass

    return possible_uniprot_ids


def fetch_structure_prediction(wdir, ref_species, gene, possible_uniprot_ids):
    """Finds and extracts structure prediction model from AlphaFold database for the species"""

    # make folder to host the structure prediction
    structure_path = wdir + "Structures/" + gene + "_" + ref_species
    if not os.path.isdir(structure_path):
        os.makedirs(structure_path)

    # remove all non pdb files and all hidden files
    clear_structure_files(structure_path)

    # find prediction model based on possible uniprot ids for the gene
    message = "\nLooking for structure in AlphaFold database for: %s ...\n" % gene
    logging.info(message)

    retrieved_uniprot_id = "No Uniprot ID retrieved"
    model_seq = ""
    gene_filepath = wdir + gene + "_" + ref_species + ".pdb"

    valid_uniprot_ids = []
    for uniprot_id in possible_uniprot_ids:
        model_version = "-F1-model_v3.pdb"  # AlphaFold database updated to v3 since writing this module
        url = "https://alphafold.ebi.ac.uk/files/AF-" + uniprot_id + model_version
        try:
            r = requests.get(url)
            r.raise_for_status()  # Raises an HTTPError if the return code is 4xx or 5xx
            with open(gene_filepath, "wb") as f:
                f.write(r.content)
            retrieved_uniprot_id = uniprot_id
            valid_uniprot_ids.append(retrieved_uniprot_id)
        except HTTPError:
            pass

    # more than one uniprot id is associated with Alpha Fold address
    # e.g. Spc25 and Spcs2 (alt. name Spc25)
    if len(valid_uniprot_ids) > 1 :
        message = "...WARNING... : More than one valid structure prediction %s detected for: %s\n" \
                  % (gene, valid_uniprot_ids)
        logging.info(message)

    if not os.path.isfile(gene_filepath):
        message = "...WARNING... : Structure prediction for gene: %s could not be fetched " \
                  "-> Cannot overlay FREEDA results onto a 3D structure\n" % gene
        logging.info(message)
        with open(structure_path + "/model_incompatible.txt", "w") as f:
            f.write("Either no reviewed entry is present in searched Ensembl release for the gene name provided or "
                    "no model has been found in AlphaFold database. Cannot overlay FREEDA results onto a 3D structure.")
        uniprot_id = None

        return model_seq, uniprot_id

    else:
        # get model sequence and length
        with open(gene_filepath, 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                model_seq = record.seq

        message = "Based on uniprot id: %s structure prediction for gene: %s has been found (%s aa)" \
                  % (retrieved_uniprot_id, gene, len(model_seq))
        logging.info(message)
        shutil.move(gene_filepath, structure_path + "/" + gene + "_" + ref_species + ".pdb")

        return model_seq, retrieved_uniprot_id


def validate_gene_names(ref_species, all_genes, all_genes_ensembl, ensembl):
    """Checks if user provided valid gene names -> find in ensembl object"""

    all_names_valid = True

    absent_names = []

    # modified 03/09/2023
    for gene in all_genes:
        if ref_species == "Dme":
            # make sure its a valid gene name
            try:
                gene_biotype = ensembl.gene_by_id(gene).biotype
                gene_name = ensembl.gene_by_id(gene).gene_name

                if gene_biotype != "protein_coding":
                    message = "... FATAL_ERROR... : %s is %s instead of a protein coding gene\n" \
                              "    -> check gene name here: ensembl.org \n" \
                              "        -> exiting the pipeline now...\n" % (gene, ensembl.gene_by_id(gene).biotype)
                    logging.info(message)
                    protein_coding = False
                    return protein_coding

                elif len(gene_name) == 1:
                    message = "...FATAL_ERROR... : %s encodes a gene called '%s' which is a common name\n" \
                              "    -> cannot reliably fetch input data \n" \
                              "        -> exiting the pipeline now...\n" % (gene, ensembl.gene_by_id(gene).gene_name)
                    logging.info(message)
                    unique_name = False
                    return unique_name

                elif gene_name not in all_genes_ensembl:
                    all_names_valid = False
                    absent_names.append(gene)

            except ValueError:
                message = "... FATAL_ERROR... : %s is NOT a valid gene name\n" \
                          "   -> check gene name here: ensembl.org \n" \
                          "        -> exiting the pipeline now...\n" % gene
                logging.info(message)
                protein_coding = False
                return protein_coding

        if ref_species != "Dme" and gene not in all_genes_ensembl:
            all_names_valid = False
            absent_names.append(gene)

    if not all_names_valid:
        message = "... FATAL_ERROR... : Gene names %s do not exist in reference assembly\n" \
                  "    -> check gene name here: ensembl.org \n" \
                  "        -> exiting the pipeline now...\n" % absent_names
        logging.info(message)

    return all_names_valid


def generate_ref_genome_object(wdir, ref_species):
    """Generates a reference Genome object using pyensembl as a wrapper for ensembl database"""

    # default is mouse for now
    if ref_species == "Mm":
        ref_genome_name = "MusMusculus_genome"
        species = "mus musculus"
        release = 104
        ref_genome_contigs_dict = genomes_preprocessing.get_ref_genome_contigs_dict(ref_species)

    elif ref_species == "Rn":
        ref_genome_name = "RattusNorvegicus_genome"
        species = "rattus norvegicus"
        release = 104
        ref_genome_contigs_dict = genomes_preprocessing.get_ref_genome_contigs_dict(ref_species)

    elif ref_species == "Hs":
        ref_genome_name = "HomoSapiens_genome"
        species = "homo sapiens"
        release = 104
        ref_genome_contigs_dict = genomes_preprocessing.get_ref_genome_contigs_dict(ref_species)

    elif ref_species == "Fc":
        ref_genome_name = "FelisCatus_genome"
        species = "felis catus"
        release = 90
        ref_genome_contigs_dict = genomes_preprocessing.get_ref_genome_contigs_dict(ref_species)

    elif ref_species == "Cf":
        ref_genome_name = "CanisFamiliaris_genome"
        species = "canis familiaris"
        release = 99   # changed from 90 on 03/21/2023
        ref_genome_contigs_dict = genomes_preprocessing.get_ref_genome_contigs_dict(ref_species)

    elif ref_species == "Gg":
        ref_genome_name = "GallusGallus_genome"
        species = "gallus gallus"
        release = 94
        ref_genome_contigs_dict = genomes_preprocessing.get_ref_genome_contigs_dict(ref_species)

    elif ref_species == "Dme":
        ref_genome_name = "DrosophilaMelanogaster_genome"
        species = "fruit fly"  # temporarily to avoid crashing on old pyensemnbl release
        release = 107
        ref_genome_contigs_dict = genomes_preprocessing.get_ref_genome_contigs_dict(ref_species)

    # make sure ref species genome (reference genome) is present
    ref_genomes_path = wdir + "Reference_genomes/"
    # check if reference genome is present -> exit if not
    ref_genome_present = tblastn.check_genome_present(wdir, ref_species, ref_genomes_path,
                                                      ref_genome_name, ref_genome=True)

    logging.getLogger("pyensembl").setLevel(logging.WARNING)  # disables logging from pyensembl
    # get ref assembly database
    ensembl = pyensembl.EnsemblRelease(release, species)
    ensembl.download()  # this is suppose to bypass installing the release from outside python
    ensembl.index()  # this is suppose to bypass installing the release from outside python

    # define biotype (FREEDA deals only with protein coding sequences)
    biotype = "Protein coding"

    # get all gene names available (list)
    all_genes_ensembl = ensembl.gene_names()
    all_gene_ids_ensembl = ensembl.gene_ids()

    return ref_genome_present, ensembl, ref_species, ref_genomes_path, ref_genome_contigs_dict, \
           biotype, all_genes_ensembl, all_gene_ids_ensembl


def extract_input(wdir, ref_species, ref_genomes_path, ref_genome_contigs_dict,
                  ensembl, biotype, gene, model_seq, uniprot_id):
    """Extracts all input sequences required by FREEDA from the indicated reference genome (ref_species)"""

    if ref_species == "Mm":
        ref_genome_name = "MusMusculus_genome"
    elif ref_species == "Rn":
        ref_genome_name = "RattusNorvegicus_genome"
    elif ref_species == "Hs":
        ref_genome_name = "HomoSapiens_genome"
    elif ref_species == "Fc":
        ref_genome_name = "FelisCatus_genome"
    elif ref_species == "Cf":
        ref_genome_name = "CanisFamiliaris_genome"
    elif ref_species == "Gg":
        ref_genome_name = "GallusGallus_genome"
    elif ref_species == "Dme":
        ref_genome_name = "DrosophilaMelanogaster_genome"

    input_correct = False
    model_matches_input = False
    microexon_present = False
    microexons = []

    blast_input_path = wdir + "Blast_input/"
    coding_sequence_input_path = wdir + "Coding_sequences/"
    gene_input_path = wdir + "Genes/"
    exons_input_path = wdir + "Exons/"

    # check if provided gene.txt list has valid gene names (if list used not GUI)
    # need to convert FlyBase ID to gene name for Drosophila
    if ref_species == "Dme":
        all_genes = {g for g in ensembl.gene_ids()}

    if ref_species != "Dme":
        all_genes = {g for g in ensembl.gene_names()}

    if gene not in all_genes:
        input_correct = False
        message = "\n...WARNING... : Reference genome lacks gene name: "'%s'"\n" % gene
        logging.info(message)
        return input_correct, model_matches_input, microexon_present, microexons

    # find coding sequence
    transcript, \
    selected_transcript_id, \
    gene_id, contig, strand, \
    UTR_5, UTR_3, \
    cds_sequence_expected, \
    matching_length = extract_cds(ensembl, ref_species, coding_sequence_input_path, gene, biotype, model_seq)

    # provided valid gene name does not code for a valid protein -> exit pipeline
    if selected_transcript_id is None:
        return input_correct, model_matches_input, microexon_present, microexons

    # find reference protein sequence
    model_matches_input = extract_protein(wdir, ref_species, blast_input_path, gene,
                                          transcript, model_seq, matching_length, uniprot_id)

    # find reference gene sequence
    extract_gene(wdir, ref_species, gene_input_path, ensembl, contig, strand, gene_id, ref_genomes_path, ref_genome_name,
                 ref_genome_contigs_dict, gene, transcript, UTR_5, UTR_3, cds_sequence_expected)

    # find reference exons sequence
    input_correct, microexon_present, microexons, missing_bp_list = extract_exons(wdir, ref_species, gene,
                                                                                  exons_input_path, ref_genomes_path,
                                                                                  ref_genome_name, contig, strand,
                                                                                  transcript, ref_genome_contigs_dict,
                                                                                  UTR_5, UTR_3, cds_sequence_expected)

    # trim reference exons and coding sequence if microexons were detected
    if microexon_present:
        input_correct = correct_for_microexons(wdir, ref_species, gene, microexons, missing_bp_list, transcript)

    # check if all expected reference exons are found in the expected reference gene sequnce
    if not check_exons_gene_compatibility(wdir, ref_species, gene):
        input_correct = False

    # sometimes poorly curated sequences would have masked residues "N" -> cannot analyse that
    if "N" in cds_sequence_expected:
        message = "\n...WARNING... : There are masked bases in the reference coding sequence -> cannot analyse"
        logging.info(message)
        input_correct = False

    # check for repetitive regions in coding sequence
    if not check_repetitive_regions(wdir, ref_species, gene):
        input_correct = False

    return input_correct, model_matches_input, microexon_present, microexons


def trim_long_gene(wdir, ref_species, gene, UTR_5, UTR_3):
    """Maps UTRs within the gene sequence of the reference gene; eliminates if > 2000bp"""

    path_to_gene = wdir + "Genes/" + gene + "_" + ref_species + "_gene.fasta"

    # check if UTRs are not too short (too many off target matches)

    # get gene sequence
    with open(path_to_gene, "r") as f:
        file = f.readlines()
        header = file[0]
        seq = file[1].upper()

    # short UTRs may align randomly
    if len(UTR_5) > 50:
        # use only the first 50 bp to map UTR
        UTR_5 = UTR_5.upper()[0:51]
        # map 5' UTR using sliding window
        for position in range(len(seq)-len(UTR_5)+1):
            if seq[position:position+len(UTR_5)] == UTR_5:
                # check if there is room to trim
                if position > 2000:
                    message = "\n...NOTE... : Reference locus for gene : %s was trimmed at 5' (performance)" % gene
                    logging.info(message)
                    trim_position = position - 2000
                    print("Length trimmed at 5' : %s" % str(trim_position))
                    # trim the gene sequence at its 5'
                    seq = seq[trim_position::]
                    break
                break

    # short UTRs may align randomly
    if len(UTR_3) > 50:
        # use only the last 50 bp to map UTR
        UTR_3 = UTR_3.upper()[len(UTR_3)-50::]
        # map 3' UTR using sliding window
        locus_length_3 = len(seq)
        for position in range(len(seq)-len(UTR_3)+1):
            # represents how much sequence is left
            locus_length_3 -= 1
            if seq[position:position+len(UTR_3)] == UTR_3:
                # check if there is room to trim
                if locus_length_3 - position - len(UTR_3) > 2000:
                    message = "\n...NOTE... : Reference locus for gene : %s was trimmed at 3' (performance)" % gene
                    logging.info(message)
                    trim_position = position + len(UTR_3) + 2000
                    print("Length trimmed at 3' : %s" % str(len(seq)-trim_position))
                    # trim the gene sequence at its 3'
                    seq = seq[0:trim_position+1]
                    break
                break

    # overwrite the gene sequence
    with open(path_to_gene, "w") as f:
        f.write(header)
        f.write(seq)


def check_exons_gene_compatibility(wdir, ref_species, gene):
    """Checks if reference exons pooled from ensembl are present in the corresponding reference gene"""

    path_to_exons = wdir + "Exons/" + gene + "_" + ref_species + "_exons.fasta"
    path_to_gene = wdir + "Genes/" + gene + "_" + ref_species + "_gene.fasta"

    with open(path_to_gene, "r") as f:
        file = f.readlines()
        for line in file:
            if not line.startswith(">") and not line.startswith("\n"):
                gene_seq = line.rstrip("\n")
                break

    with open(path_to_exons, "r") as f:
        file = f.readlines()
        for line in file:
            if not line.startswith(">") and not line.startswith("\n"):
                if line.rstrip("\n") not in gene_seq:
                    message = "\n...FATAL ERROR... : At least one expected ref exon ABSENT " \
                              "in expected ref gene sequence"
                    logging.info(message)
                    return False

        return True


def extract_protein(wdir, ref_species, blast_input_path, gene, transcript, model_seq, matching_length, uniprot_id):
    """Extracts protein sequence based on Transcript object and compares with structure model sequence"""

    model_matches_input = False
    structure_path = wdir + "Structures/" + gene + "_" + ref_species

    # find protein sequence
    protein_sequence = transcript.protein_sequence
    # get and save the sequence
    parse_sequence(ref_species, blast_input_path, protein_sequence,
                   gene, transcript, strand=None, sequence_type="protein")
    # check if protein length and sequence match that of the model
    if matching_length is True and model_seq == protein_sequence:
        model_matches_input = True
        with open(structure_path + "/model_matches_input_seq.txt", "w") as f:
            f.write("Model is based on an identical protein sequence as blast input.")
            f.write("\nUniprot ID : %s" % uniprot_id)

    else:
        with open(structure_path + "/model_incompatible.txt", "w") as f:
            f.write("Model sequence does not match the protein sequence used for blast input. "
                    "Cannot overlay FREEDA results onto a 3D structure.")
            f.write("\nUniprot ID : %s" % uniprot_id)

    return model_matches_input


def extract_exons(wdir, ref_species, gene, exons_input_path, ref_genomes_path, ref_genome_name,
                  contig, strand, transcript, ref_genome_contigs_dict, UTR_5, UTR_3, cds_sequence_expected):
    """Extracts exoms sequence based on Transcript object"""

    # get start and stop for each exon
    exon_intervals = transcript.exon_intervals

    start_codon_offset = transcript.first_start_codon_spliced_offset

    # get sequence of each exon
    exons_filename = gene + "_" + ref_species + "_exons.fasta"
    exons_file = open(exons_input_path + exons_filename, "w")

    number = 0
    header = ">"
    coding_exons = {}
    UTR_5_length = len(UTR_5)
    UTR_3_length = len(UTR_3)
    removed_non_coding_from_start = False
    removed_non_coding_from_end = False
    start_codon_present = False
    stop_codon_present = False
    nr_of_removed_non_coding_exons_start = 0
    nr_of_removed_non_coding_exons_end = 0
    cds_from_exons = ""
    microexon_present = False
    input_correct = True
    microexons = []
    missing_bp_list = []

    # make a true copy of exon coordinates
    exon_intervals_copy = exon_intervals[::]

    # Get first coding exon
    for exon in exon_intervals_copy:

        if exon[1]-exon[0] > UTR_5_length:
            break

        # trim off the non-coding exons containing UTR_5
        if exon[1]-exon[0] <= UTR_5_length:
            removed_non_coding_from_start = True
            nr_of_removed_non_coding_exons_start += 1
            exon_intervals = exon_intervals[1:]
            trimmed_UTR_5_length = UTR_5_length - (exon[1]-exon[0])
            start_codon_offset = trimmed_UTR_5_length
            UTR_5_length = trimmed_UTR_5_length

    # Get last coding exon
    for exon in reversed(exon_intervals_copy):

        if exon[1]-exon[0] > UTR_3_length:
            break

        # trim off the non-coding exons containing UTR_3
        if exon[1]-exon[0] <= UTR_3_length:
            removed_non_coding_from_end = True
            nr_of_removed_non_coding_exons_end += 1
            exon_intervals = exon_intervals[:-1]
            trimmed_UTR_3_length = UTR_3_length - (exon[1]-exon[0])
            UTR_3_length = trimmed_UTR_3_length

    # Get all coding exons
    for exon in exon_intervals:
        number += 1
        start = exon[0] - 1
        stop = exon[1]

        if number == 1 and strand == "+":
            if removed_non_coding_from_start is False:
                start = start + start_codon_offset
            else:
                start = start + start_codon_offset - nr_of_removed_non_coding_exons_start

        if number == 1 and strand == "-":
            if removed_non_coding_from_start is False:
                stop = stop - start_codon_offset
            else:
                stop = stop - start_codon_offset + nr_of_removed_non_coding_exons_start

        if number == len(exon_intervals) and strand == "+":
            if removed_non_coding_from_end is False:
                stop = stop - UTR_3_length
            else:
                stop = stop - UTR_3_length + nr_of_removed_non_coding_exons_end

        if number == len(exon_intervals) and strand == "-":
            if removed_non_coding_from_end is False:
                start = start + UTR_3_length
            else:
                # this does not behave well when more than 1 non-coding exon at 3' e.g. FBgn0016976
                # ADDED "nr_of_removed_non_coding_exons_end" instead of "-1" on 03/09/2023
                start = start + UTR_3_length - nr_of_removed_non_coding_exons_end

        exon_fasta_sequence = get_single_exon(ref_species, gene, ref_genome_contigs_dict,
                ref_genomes_path, ref_genome_name, rules, contig, strand, number, start, stop)
        coding_exons[number] = exon_fasta_sequence

    for nr, exon_sequence in coding_exons.items():

        if len(exon_sequence) < 18:
            with open(wdir + "Structures/" + gene + "_" + ref_species + "/microexons.txt", "a") as f:
                f.write(str(nr) + "\n")
            microexon_present = True
            missing_bp_list.append(len(exon_sequence))
            microexons.append((nr, str(len(exon_sequence)) + "bp"))
            message = "\n...WARNING...: Exon %s in: %s is a microexon (%sbp; hard to align) " \
                    "-> eliminated from exons and cds" % (str(nr), gene, len(exon_sequence))
            logging.info(message)
            # check if the last exon (microexon) had a start codon before its trimmed (ex. Pot1a)
            if nr == 1 and exon_sequence.startswith("ATG"):
                start_codon_present = True
            # check if the last exon (microexon) had a stop codon before its trimmed (ex. Numa1)
            if nr == len(coding_exons) and exon_sequence[-3:] in ["TGA", "TAG", "TAA"]:
                stop_codon_present = True
            continue

        else:
            cds_from_exons += exon_sequence

            if nr == 1 and exon_sequence.startswith("ATG"):
                start_codon_present = True
            if nr == len(coding_exons) and exon_sequence[-3:] in ["TGA", "TAG", "TAA"]:
                stop_codon_present = True

        header = ">" + transcript.name + "_" + ref_species + "_exon_" + str(nr)
        exons_file.write(header + "\n")
        exons_file.write(exon_sequence + "\n")

    if start_codon_present is False:
        message = "\nFirst exon in: %s is missing a START codon!!!\n" % header.split("_")[0].replace(">","")
        logging.info(message)
        input_correct = False

    if stop_codon_present is False:
        message = "\nLast exon in: %s is missing a STOP codon!!!\n" % header.split("_")[0].replace(">","")
        logging.info(message)
        input_correct = False

    if microexon_present is False and cds_from_exons != cds_sequence_expected:
        message = "\nExons FAILED to assemble expected CDS for: %s\n" % header.split("_")[0].replace(">","")
        logging.info(message)
        input_correct = False

    exons_file.close()

    return input_correct, microexon_present, microexons, missing_bp_list


def get_single_exon(ref_species, gene, ref_genome_contigs_dict, ref_genomes_path,
                    ref_genome_name, rules, contig, strand, number, start, stop):
    """Returns a header and a string representation of a current exon."""

    header = gene + "_" + ref_species + "_exon_" + str(number)
    exon_bed_filename = header + ".bed"
    exon_fasta_filename = header + ".fasta"
    with open(exon_bed_filename, "w") as b:
        b.write(ref_genome_contigs_dict[contig] + "\t")
        b.write(str(start) + "\t")
        b.write(str(stop))

    bed_object = pybedtools.BedTool(exon_bed_filename)
    bed_object = bed_object.sequence(fi = ref_genomes_path + ref_genome_name + ".fasta")
    # makes a file and writes a fasta seq into it
    exon_fasta_sequence = bed_object.save_seqs(exon_fasta_filename)

    with open(exon_fasta_filename, "r") as f:
        file = f.readlines()
        exon_fasta_sequence = file[1].rstrip("\n").upper()

    if strand == "-":
        complement = ""
        for base in exon_fasta_sequence:
            complement = rules[base] + complement
        exon_fasta_sequence = complement

    os.remove(exon_bed_filename)
    os.remove(exon_fasta_filename)

    return exon_fasta_sequence


def extract_gene(wdir, ref_species, gene_input_path, ensembl, contig, strand, gene_id,
                 ref_genomes_path, ref_genome_name, ref_genome_contigs_dict, gene,
                 transcript, UTR_5, UTR_3, cds_sequence_expected):
    """Extracts gene sequence based on Genome object"""

    # list all genes in contig
    genes = ensembl.genes(contig)

    # create a Gene object
    if ref_species == "Dme":
        gene_obj = [g for g in genes if g.gene_id == gene_id][0]
    if ref_species != "Dme":
        gene_obj = [g for g in genes if g.gene_id == gene_id[0]][0]

    # get gene name
    gene_name = gene_obj.gene_name
    # get gene starting position
    start = gene_obj.start - 200   # extend flanking regions missing in some genes (e.g. TLR5 Gg)
    # get gene end position
    end = gene_obj.end + 200  # extend flanking regions missing in some genes (e.g. TLR5 Gg)
    # get gene sequence

    # make a bed and fasta file for gene (add underscore to differenciate from other handles)
    gene_bed_filename = "_" + gene_name + ".bed"

    matched_contig = ref_genome_contigs_dict[str(contig)]

    # make a bedtool file
    with open(gene_bed_filename, "w") as b:
        b.write(str(matched_contig) + "\t" + str(start) + "\t" + str(end))

    # parse that bedtool file using BedTools object
    bed_object = pybedtools.BedTool(gene_bed_filename)
    gene_fasta_sequence = bed_object.sequence(fi=ref_genomes_path + ref_genome_name + ".fasta")

    os.remove(gene_bed_filename)

    parse_sequence(ref_species, gene_input_path, gene_fasta_sequence,
                   gene, transcript, strand, sequence_type="gene")

    # last bp is often missing in gene where UTR_3 is missing -> last exon will not be called unless fixed
    if len(UTR_3) == 0:

        if ref_species == "Dme":
            gene_name = gene_id

        with open(wdir + "Genes/" + gene_name + "_" + ref_species + "_gene.fasta", "r") as f:
            file = f.readlines()

            if file[1].endswith("TG") or file[1].endswith("TA"):  # check if gene sequence ends with a partial STOP
                message = "\n...WARNING... : 3' UTR not detected in %s and gene ends with a partial STOP " \
                      "or exactly at STOP -> added missing bp" % gene
                logging.info(message)

                with open(wdir + "Genes/" + gene_name + "_" + ref_species + "_gene.fasta", "w") as w:
                    w.write(file[0].rstrip("\n"))  # added on 11/13/2021

                    # add the missing bp based on STOP from cds AND ADD A PLACEHOLDER BP TO FACILITATE EXON CALLING
                    if cds_sequence_expected[-3:] == "TGA" or cds_sequence_expected[-3:] == "TAA":
                        w.write("\n" + file[1] + "A" + "A")  # TGA or TAA + placeholder bp
                    else:
                        w.write("\n" + file[1] + "G" + "G")  # TAG + placeholder bp

    trim_long_gene(wdir, ref_species, gene, UTR_5, UTR_3)


def parse_sequence(ref_species, output_path, fasta_sequence, gene, transcript, strand, sequence_type):
    """Parses either a pybedtools.bedtool.BedTool object or string into a file"""

    sequence = ""
    header = ""

    # check if its a bedtool object
    if isinstance(fasta_sequence, pybedtools.bedtool.BedTool):
        file = open(fasta_sequence.seqfn)
        content = file.readlines()
        header = ">" + transcript.name + "_" + ref_species + "_" + sequence_type
        sequence = content[1].upper().rstrip("\n")
        file.close()

    # or its a plain string
    if not isinstance(fasta_sequence, pybedtools.bedtool.BedTool):
        header = ">" + transcript.name + "_" + ref_species + "_" + sequence_type
        sequence = fasta_sequence

    # or plain string on reverse string (gene only)
    if sequence_type == "gene" and strand == "-":
        complement = ""
        for base in sequence:
            complement = rules[base] + complement
        sequence = complement

    # write it into a file
    filename = gene + "_" + ref_species + "_" + sequence_type
    with open(output_path + filename + ".fasta", "w") as f:
        f.write(header + "\n")
        f.write(sequence)


def extract_cds(ensembl, ref_species, coding_sequence_input_path, gene, biotype, model_seq):
    """Extracts coding sequence by creating Transcript object"""

    # ADDED 03/09/2023 to prevent syntax errors with various symbols
    if ref_species == "Dme":
        all_transcripts_ids = ensembl.transcript_ids_of_gene_id(gene)
    if ref_species != "Dme":
        all_transcripts_ids = ensembl.transcript_ids_of_gene_name(gene)

    all_transcripts_dict = {}
    selected_transcript_id = None
    transcript = None
    UTR_5 = None
    UTR_3 = None
    gene_id = None
    contig = None
    strand = None
    matching_length = False
    cds_sequence_expected = ""

    # find all possible transcripts
    for t in all_transcripts_ids:

        all_transcripts_dict[t] = []

        # find gene id
        if ref_species != "Dme":
            gene_id = ensembl.gene_ids_of_gene_name(gene)
        if ref_species == "Dme":
            gene_id = gene

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
        transcript = pyensembl.Transcript(t, transcript_name,
            contig, start, end, strand, biotype, gene_id, ensembl, support_level=None)

        # try to get cds -> if "start_codon" does not exist -> KeyError and then ValueError occur
        try:
            cds_sequence_expected = transcript.coding_sequence
        except ValueError:
            cds_sequence_expected = ""
            all_transcripts_dict[t].append(cds_sequence_expected)
        finally:
            all_transcripts_dict[t].append(cds_sequence_expected)

        # get length of cds
        if cds_sequence_expected is not None:
            length = len(cds_sequence_expected)
            all_transcripts_dict[t].append(length)

    # compare length of translated coding sequence with model amino acid sequence length
    likely_transcripts = {}  # collect likely transcripts
    for t, features in all_transcripts_dict.items():

        # some transcrtipts do not have an expected cds sequence (non-coding transcripts)
        if features[-1] is not None:

            cds_length_aa = (features[-1]-3)/3
            if len(model_seq) == cds_length_aa:
                likely_transcripts[t] = features
                matching_length = True

    # select the most likely transcript by taking the lowest ensembl number (usually the canonical transcript)

    # Drosophila transcripts dont have names followed by "-201", "-202"... etc conncention
    # but rather "-RA", "-RB" etc. so the ensembl_number cannot be an integer

    if ref_species == "Dme":

        number = "A"
        for t, features in likely_transcripts.items():
            # need to pick the last element cose some gene names have "-R" as part of their name
            ensembl_number = str(features[3].split("-R")[-1])
            if ensembl_number < number or number == "A":
                number = ensembl_number
                selected_transcript_id = t

    if ref_species != "Dme":

        number = float("inf")
        for t, features in likely_transcripts.items():
            ensembl_number = int(features[3].split("-")[1])
            if ensembl_number < number:
                number = ensembl_number
                selected_transcript_id = t

    # if there is no transcript of preference, pick the one with longest cds (not recommended)
    if matching_length is False:
        length = 0

        for t, features in all_transcripts_dict.items():

            # some transcrtipts do not have an expected cds sequence (non-coding transcripts)
            if features[-1] is not None:
                if features[-1] > length:
                    selected_transcript_id = t
                    length = features[-1]

    # some gene names in ensembl do not code proteins -> exit pipeline
    if not selected_transcript_id:
        message = "\n...FATAL ERROR... : No reliable coding sequence annotation detected for gene %s" % gene
        logging.info(message)
        return transcript, selected_transcript_id, gene_id, contig, strand, \
               UTR_5, UTR_3, cds_sequence_expected, matching_length

    # unpack features of the selected transcript
    gene_id, contig, strand, transcript_name, start, end, \
    cds_sequence_expected, length = all_transcripts_dict[selected_transcript_id]
    # rebuild Transcript object based on the selected transcript
    transcript = pyensembl.Transcript(selected_transcript_id, transcript_name,
            contig, start, end, strand, biotype, gene_id, ensembl, support_level=None)

    # check if given transcript has an annotated start codon
    start_codon_present = transcript.contains_start_codon
    # get 5'UTR sequence
    UTR_5 = transcript.five_prime_utr_sequence.upper()
    # check if given transcript has an annotated stop codon
    stop_codon_present = transcript.contains_stop_codon
    # get 3'UTR sequence
    UTR_3 = transcript.three_prime_utr_sequence.upper()
    # get and save the sequence
    parse_sequence(ref_species, coding_sequence_input_path, cds_sequence_expected,
                   gene, transcript, strand, sequence_type="cds")

    # check for START and STOP codons
    if start_codon_present and stop_codon_present:
        return transcript, selected_transcript_id, gene_id, \
               contig, strand, UTR_5, UTR_3, \
               cds_sequence_expected, matching_length
    else:
        message = "Chosen transcript for gene: %s does not have START and STOP annotated" % gene
        logging.info(message)


def clear_structure_files(structure_path):
    """Clears all files that are not pdb"""

    model_file_list = os.listdir(structure_path)
    all_files = [file for file in model_file_list]

    # remove all files
    for file in all_files:
        try:
            os.remove(structure_path + "/" + file)
        except FileNotFoundError:
            pass

    return


def check_microexons(wdir, gene, ref_species):
    """Checks if microexons were found during automatic input extraction"""

    path_to_model_info = wdir + "Structures/" + gene + "_" + ref_species
    if os.path.isfile(path_to_model_info + "/microexons.txt"):
        with open(path_to_model_info + "/microexons.txt", "r") as f:
            file = f.readlines()
            microexons = [int(exon.rstrip("\n")) for exon in file]
    else:
        microexons = []

    return microexons


def check_repetitive_regions(wdir, ref_species, gene):
    """Checks if there are large repetitive regions within the coding sequence"""

    repeat_length = 80
    path_to_cds = wdir + "Coding_sequences/" + gene + "_" + ref_species + "_cds.fasta"

    with open(path_to_cds) as f:
        seq = f.readlines()[1]

        # collect sequences using sliding window
        possible_repeats = {}
        for i in range(len(seq) - repeat_length):
            if seq[i:i + repeat_length + 1] not in possible_repeats:
                repeat = seq[i:i + repeat_length + 1]
                possible_repeats[repeat] = 1
            else:
                message = "\n...FATAL ERROR... : Repetitive coding sequence detected in %s (min 80bp repeat)" \
                          "\n   -> cannot reliably analyze this gene" % gene
                logging.info(message)
                return False

    # no repetitive sequences detected
    return True


def get_gene_names(wdir):  # THIS FUNCTION IS ONLY FOR TESTING
    """Gets sample of gene names from ensembl using pyensembl package"""

    import re
    release = 104
    species = "homo sapiens"
    ensembl = pyensembl.EnsemblRelease(release, species)
    ensembl.download()  # this is suppose to bypass installing the release from outside python
    ensembl.index()  # this is suppose to bypass installing the release from outside python

    all_genes_ensembl = ensembl.gene_names()
    with open(wdir + "proteins_temp.txt", "w") as f:

        genes = []
        for i in range(5000, len(all_genes_ensembl), 150):
            gene = all_genes_ensembl[i]
            if not gene.startswith("MIR") \
                    and not gene.startswith("LINC") \
                    and "-" not in gene \
                    and "PS" not in gene \
                    and "OS" not in gene \
                    and "IK" not in gene \
                    and "SNO" not in gene \
                    and "GM" not in gene \
                    and "RNA" not in gene \
                    and not re.search(r"P\d{0,2}$", gene):  # pseudogenes often or dont have transcripts
                genes.append(gene)
                f.write(gene + "\n")