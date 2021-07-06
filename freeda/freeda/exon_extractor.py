#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 11:00:54 2021

@author: damian
"""

from freeda import folder_generator
from freeda import name_finder
from freeda import gene_and_cds_reader
from freeda import genome_indexer
from freeda import matches_generator
from freeda import matches_processor
from freeda import msa_aligner
from freeda import msa_analyzer
import datetime
import glob
import time
import logging
import shutil


def analyse_blast_results(wdir, blast_path, original_species, t, mafft_path):   
    
    start_time = time.time()

    #if mafft_path == None:
    #    mafft_path = "/Users/damian/anaconda3/bin/mafft"
    
    day = datetime.datetime.now().strftime("-%m-%d-%Y-%H-%M")
    result_path = wdir + "Results" + day + "/"
    
    folder_generator.generate_folders(result_path, t)
    
    # initiate log file to record PAML analysis by reseting the handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    log_filename = "FREEDA" + day + ".log"
    logging.basicConfig(filename=log_filename, level=logging.INFO, format="%(message)s")

    # make a list of paths with blast tables
    all_blasts = [blast for blast in glob.glob(blast_path + "*.txt")]
    
    # get path for a single blast table while removing it from the list
    for path in all_blasts:
        match_path = path
        # generate protein and genome names and make it global
        protein_name, genome_name = name_finder.get_names(match_path)
        # find cds and gene for this path (Mm_exons and microexons NOT USED HERE)
        cds, gene, Mm_exons, microexons, expected_exons = gene_and_cds_reader.find_gene_and_cds(wdir, protein_name, original_species)
        # index given genome
        genome_index = genome_indexer.index_genome_database(wdir, genome_name)
        # generate matches dataframe
        matches = matches_generator.generate_matches(match_path, t, protein_name, genome_name, genome_index)
        # process the final dataframe
        MSA_path = matches_processor.process_matches(wdir, matches, cds, gene, t, result_path, protein_name, genome_name, genome_index)        
        # run MAFFT on all the MSA and write them into files
        msa_aligner.run_MAFFT(MSA_path, mafft_path)
        # return potential exons for a current protein in current genome
        msa_analyzer.analyse_MSA(wdir, MSA_path, protein_name, genome_name, mafft_path, result_path, Mm_exons, microexons, expected_exons)
        # mark that this blast result has been analysed
        message = "\nFinished running protein: '%s' from genome: '%s'\n" \
            % (protein_name, genome_name)
        logging.info(message)
        print(message)
        
    # generate a list of all files in the working directory
    from os import listdir
    from os.path import isfile, join
    all_files = [f for f in listdir(wdir) if isfile(join(wdir, f))]
    
    # generate a list of all protein names
    with open("proteins.txt", "r") as f:
        all_proteins = [protein.rstrip("\n") for protein in f.readlines()]
    
    # add original_species cds for a give protein
    for protein in all_proteins:
        if protein + ".fasta" in all_files:
            seq = get_original_cds(wdir, protein, original_species)
            header = ">" + protein + "_" + original_species
            with open(protein + ".fasta", "r+") as file:
                content = file.read()
                file.seek(0, 0)
                file.write(header.rstrip('\r\n') + '\n' + seq + '\n' + content)
    
    
    # analyse_final_cds(original_species, result_path, mafft_path) # SHOULD BE RUN SEPARATELY
    
    # mark the end of the analysis
    message = ("Analysis completed in %s minutes or %s hours" % \
               ((time.time() - start_time)/60, 
                (time.time() - start_time)/60/60))
    print(message)
    logging.info(message)
    
    # move the log file into result folder
    shutil.move(log_filename, result_path)
    
    return result_path


def get_original_cds(wdir, protein_name, original_species):

    # open according cds fasta file
    with open(wdir + "Coding_sequences/" + protein_name + "_" + original_species + "_cds.fasta", "r") as f:
        sequence = ""
        cds = f.readlines()
        for line in cds[1:]:
            sequence = sequence + line.rstrip("\n")
    return sequence