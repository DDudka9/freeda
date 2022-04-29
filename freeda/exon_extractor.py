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


def analyze_blast_results(wdir, blast_output_path, ref_species, t, all_genes, all_genomes,
                          final_excluded_species=None, gui=None, logging_window=None, all_genes_dict=None):
    """ Finds and clones exons based on blast results"""

    start_time = time.time()

    day = datetime.datetime.now().strftime("-%m-%d-%Y-%H-%M")
    result_path = wdir + "Results" + day + "/Raw_data/"
    log_filename = "FREEDA" + day + ".log"

    # does not generate folders for excluded species
    folder_generator.generate_contig_folders(result_path, all_genes, all_genomes)

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
    all_blasts = sorted(glob.glob(blast_output_path + "*.txt"), key=os.path.getmtime, reverse=False)

    # remove blast output files for excluded species
    for blast in all_blasts:
        for species, genome in final_excluded_species.items():
            if genome in blast:  # find genome name e.e. "MusSpicilegus" in absolute path of blast file
                all_blasts.remove(blast)

    # remove previous fasta files for these genes (unfinished runs)
    for gene in all_genes:
        if os.path.isfile(wdir + gene + ".fasta"):
            os.remove(wdir + gene + ".fasta")

    # get path for a single blast table while removing it from the list
    for gene in all_genes:
        for path in all_blasts:
            if gene == path.split("/")[-1].split("_")[0]:  # to avoid gene names hidden in path namespace

                # copy to the Result folder
                shutil.copy(path, result_path.replace("Raw_data/", "Results/") + "Blast_output/" + path.split("/")[-1])

                match_path = path
                # generate gene and genome names
                gene, genome_name = get_names(match_path)
                # find cds and gene for this path (ref_exons NOT USED HERE)
                cds_seq, gene_seq, ref_exons, expected_exons = fasta_reader.find_gene_and_cds(wdir,
                                                                                              gene, ref_species)
                # index given genome
                genome_index = genome_indexer.index_genome_database(wdir, genome_name)
                # generate matches dataframe
                matches = matches_generator.generate_matches(match_path, t, gene, genome_name)
                # process the final dataframe
                MSA_path = matches_processor.process_matches(ref_species, wdir, matches, cds_seq, gene_seq, result_path,
                                                             gene, genome_name, genome_index)
                # run MSA and write them into files
                msa_aligner.run_msa(MSA_path)
                # return potential exons for a current gene in current genome
                msa_analyzer.analyze_MSA(wdir, ref_species, MSA_path, gene,
                                         genome_name, ref_exons, expected_exons, all_genes_dict)
                # mark that this blast result has been analyzed
                message = "\nFinished running gene: '%s' from genome: '%s'\n" \
                    % (gene, genome_name)
                logging.info(message)
                print(message)

    # generate a list of all files in the working directory
    all_files = [f for f in os.listdir(wdir) if os.path.isfile(os.path.join(wdir, f))]

    # add ref_species cds for a give gene
    for gene in all_genes:
        if gene + ".fasta" in all_files:
            seq = get_ref_cds(wdir, gene, ref_species)
            header = ">" + gene + "_" + ref_species
            with open(gene + ".fasta", "r+") as file:
                content = file.read()
                file.seek(0, 0)
                file.write(header.rstrip('\r\n') + '\n' + seq + '\n' + content)
        shutil.move(gene + ".fasta", result_path)

    # mark the end of the analysis
    message = ("Analysis completed in %s minutes or %s hours" %
               ((time.time() - start_time)/60,
                (time.time() - start_time)/60/60))
    print(message)
    logging.info(message)

    # move the log file into result folder
    shutil.move(log_filename, result_path)

    return result_path


def check_blast_output(blast_output_path, t, all_genes):
    """Checks if at least one blast output matches for a given gene passes the blast threshold picked by the user"""

    blast_output_correct = True

    for gene in all_genes:

        blast_output_files = [file for file in os.listdir(blast_output_path) if file.startswith(gene)]
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
                          "for gene name: %s -> exiting the pipeline now..." % (t, gene))
                    return blast_output_correct

    return blast_output_correct


def get_ref_cds(wdir, gene, ref_species):
    """Reads cds for the gene in the reference species (from "Coding_sequences" folder"""

    # open according cds fasta file
    with open(wdir + "Coding_sequences/" + gene + "_" + ref_species + "_cds.fasta", "r") as f:
        sequence = ""
        cds = f.readlines()
        for line in cds[1:]:
            sequence = sequence + line.rstrip("\n")

    return sequence


def get_names(match_path):
    """Finds gene names and genome names in the matches path"""

    # isolate flag names for gene and genome from blast result filename:
    # get blast result name
    path_split = re.split(r"/", match_path)[-1]
    # get gene name
    gene = re.split(r"_", path_split)[0]
    # get genome name (in 3 steps)
    genome_file_name = re.split(r"_", path_split)[1:3]
    genome_suffix = re.split(r"\.", genome_file_name[1])[0]
    genome_name = genome_file_name[0] + "_" + genome_suffix

    return gene, genome_name

