#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 11:00:54 2021

@author: damian
"""

from freeda import folder_generator
from freeda import fasta_reader
from freeda import genome_indexer
from freeda import matches_generator
from freeda import matches_processor
from freeda import msa_aligner
from freeda import msa_analyzer
from freeda import TextHandler
import datetime
import glob
import time
import logging
import shutil
import os
import re


def analyse_blast_results(wdir, blast_output_path, ref_species, t, all_proteins, all_genomes, aligner, gui=None,
                          logging_window=None, all_proteins_dict=None):
    """ Finds and clones exons based on blast results"""

    start_time = time.time()

    day = datetime.datetime.now().strftime("-%m-%d-%Y-%H-%M")
    result_path = wdir + "Results" + day + "/"
    log_filename = "FREEDA" + day + ".log"

    folder_generator.generate_folders(result_path, all_proteins, all_genomes)

    # initiate log file to record FREEDA analysis by reseting the handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    if gui:
        # make a new handler
        text_handler = TextHandler.TextHandler(logging_window)
        # configure the new logger
        logging.basicConfig(filename=log_filename, level=logging.INFO, format="%(message)s")
        logger = logging.getLogger()
        logger.addHandler(text_handler)

    else:
        # configure the logger
        logging.basicConfig(filename=log_filename, level=logging.INFO, format="%(message)s")

    # make a list of paths with blast tables
    all_blasts = [blast for blast in glob.glob(blast_output_path + "*.txt")]

    # remove previous fasta files for these proteins (unfinished runs)
    for protein in all_proteins:
        if os.path.isfile(wdir + protein + ".fasta"):
            os.remove(wdir + protein + ".fasta")

    # get path for a single blast table while removing it from the list
    for protein in all_proteins:
        for path in all_blasts:
            if protein in path:
                match_path = path
                # generate protein and genome names and make it global
                protein_name, genome_name = get_names(match_path)
                # find cds and gene for this path (ref_exons NOT USED HERE)
                cds, gene, ref_exons, expected_exons = fasta_reader.find_gene_and_cds(wdir, protein_name, ref_species)
                # index given genome
                genome_index = genome_indexer.index_genome_database(wdir, genome_name)
                # generate matches dataframe
                matches = matches_generator.generate_matches(match_path, t, protein_name, genome_name, genome_index)
                # process the final dataframe
                MSA_path = matches_processor.process_matches(wdir, matches, cds, gene, result_path,
                                                             protein_name, genome_name, genome_index)
                # run MSA and write them into files
                msa_aligner.run_msa(MSA_path, aligner)
                # return potential exons for a current protein in current genome
                msa_analyzer.analyse_MSA(wdir, ref_species, MSA_path, protein_name,
                                         genome_name, ref_exons, expected_exons, aligner, all_proteins_dict)
                # mark that this blast result has been analysed
                message = "\nFinished running protein: '%s' from genome: '%s'\n" \
                    % (protein_name, genome_name)
                logging.info(message)
                print(message)

    # generate a list of all files in the working directory
    all_files = [f for f in os.listdir(wdir) if os.path.isfile(os.path.join(wdir, f))]

    # add ref_species cds for a give protein
    for protein in all_proteins:
        if protein + ".fasta" in all_files:
            seq = get_ref_cds(wdir, protein, ref_species)
            header = ">" + protein + "_" + ref_species
            with open(protein + ".fasta", "r+") as file:
                content = file.read()
                file.seek(0, 0)
                file.write(header.rstrip('\r\n') + '\n' + seq + '\n' + content)
        shutil.move(protein + ".fasta", result_path)

    # mark the end of the analysis
    message = ("Analysis completed in %s minutes or %s hours" %
               ((time.time() - start_time)/60,
                (time.time() - start_time)/60/60))
    print(message)
    logging.info(message)

    # move the log file into result folder
    shutil.move(log_filename, result_path)

    return result_path


def check_blast_output(blast_output_path, t, all_proteins):
    """Checks if at least one blast output matches for a given protein passes the blast threshold picked by the user"""

    blast_output_correct = True

    for gene_name in all_proteins:

        blast_output_files = [file for file in os.listdir(blast_output_path) if file.startswith(gene_name)]
        # make sure there are no hidden files
        genome_names = [file for file in blast_output_files if not file.startswith(".")]
        no_matches_above_t = 0

        for genome_name in genome_names:

            with open(blast_output_path + genome_name, "r") as f:
                file = f.readlines()

                if not [match for match in file if float(match.split("\t")[9]) > t]:
                    no_matches_above_t += 1
                    os.remove(blast_output_path + genome_name)
                    message = "\n...WARNING... : No matches above threshold : %s " \
                              "found in blast output file : %s" % (t, genome_name)
                    print(message)
                    logging.info(message)

                if no_matches_above_t > 3:
                    blast_output_correct = False
                    print("\n...FATAL ERROR... : At least 3 blast output files contain no matches above threshold : %s "
                          "for gene name: %s -> exiting the pipeline now..." % (t, gene_name))
                    return blast_output_correct

    return blast_output_correct


def get_ref_cds(wdir, protein_name, ref_species):
    """Reads cds for the protein in the reference species (from "Coding_sequences" folder"""

    # open according cds fasta file
    with open(wdir + "Coding_sequences/" + protein_name + "_" + ref_species + "_cds.fasta", "r") as f:
        sequence = ""
        cds = f.readlines()
        for line in cds[1:]:
            sequence = sequence + line.rstrip("\n")
    return sequence


def get_names(match_path):
    # isolate flag names for protein and genome from blast result filename:
    # get blast result name
    path_split = re.split(r"/", match_path)[-1]
    # get protein name
    protein_name = re.split(r"_", path_split)[0]
    # get genome name (in 3 steps)
    genome_file_name = re.split(r"_", path_split)[1:3]
    genome_suffix = re.split(r"\.", genome_file_name[1])[0]
    genome_name = genome_file_name[0] + "_" + genome_suffix
    return protein_name, genome_name