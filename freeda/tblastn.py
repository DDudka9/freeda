#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:21:53 2021

@author: damian

Runs blast using NCBI tblastn.

Fasta genome names need to be stored in "genomes.txt" file one per line using one underscore:
    xxxxx_xxxxx.fasta
    yyyyy_yyyyy.fasta
    
Gene names need to be stored in "genes.txt" file one per line:
    aaaaa
    bbbbb

Outputs tabulated text files per gene per genomes.

"""

from freeda import genomes_preprocessing, pyinstaller_compatibility
import subprocess
import os
import glob
import time
import logging


def run_blast(wdir, ref_species, all_genes, final_excluded_species=None):
    """Runs tblastn based on NCBI makedatabase routine."""

    database_path = wdir + "Genomes/"
    query_path = wdir + "Blast_input/"  # Coding_sequences/
    output_path = wdir + "Blast_output/"
    task = "blastn"  # this is to delete
    evalue = "0.01"
    form = "6 qseqid means sseqid means qstart means qend means sstart means send means evalue means " \
           "bitscore length means pident means mismatch means gapopen means qlen means slen means"

    all_genomes = genomes_preprocessing.get_names(wdir, ref_species, final_excluded_species)
    genomes = [names[1] for names in all_genomes]

    # clear Blast_output folder
    all_old_blast_output_files = glob.glob(os.path.join(output_path, "*.txt"))
    for file in all_old_blast_output_files:
        os.remove(file)

    # make sure database for each genome already exists or is successfully built
    for genome in genomes:
        genome_file_database = check_genome_present(wdir, ref_species, database_path, genome, ref_genome=False)
        # failed to build database
        if genome_file_database is False:
            return None

    # PYINSTALLER: Add path to os path variable.
    tblastn_path = pyinstaller_compatibility.resource_path("tblastn")  # from blastn

    # perform blast
    for genome in genomes:
        for gene in all_genes:
            database = database_path + genome + ".fasta"
            query = query_path + gene + "_" + ref_species + "_protein.fasta"  # from _cds.fasta
            output = output_path + gene + "_" + genome + ".txt"
            to_blast = [tblastn_path, "-db", database, "-query", query, "-evalue", evalue,  # added evalue 05_12_2022
                        "-out", output, "-outfmt", form, "-num_threads", "8"]
            #to_blast = [tblastn_path, "-db", database, "-query", query, "-task", task, "-evalue", evalue,
            #            "-out", output, "-outfmt", form, "-num_threads", "8"]
            message = "\nPerforming tblastn for gene: %s from genome: %s\n" % (gene, genome)
            logging.info(message)
            subprocess.call(to_blast)

    print("\ntblastn txt files have been generated.")
    
    return


def check_genome_present(wdir, ref_species, database_path, genome, ref_genome=False):
    """Checks if a given genome is present. Unpacks and unzips genomes downloaded from NCBI Assembly.
    Non-ncbi assemblies must be prepared as ".fasta" files conform with "genomes.txt" names.
    It also looks for reference genome if key-only argument reference_genome is invoked."""

    genome_file_database = True

    # if function called by input extractor module not tblastn module
    if ref_genome is True:
        ref_genome_path = database_path
    
    # define expected files
    zip_file = genome + ".zip"
    expected_genome_file = genome + ".fasta"
    
    # get info on all files available
    all_files = []
    for root, dirs, files in os.walk(database_path, topdown=False):
        for f in files:
            all_files.append(f)

    # there are no genomes in the folder
    if not all_files:
        all_files.append("empty")
    
    # check if genome database is present
    for file in all_files:
        # look for indices
        if expected_genome_file + ".pal" not in all_files and expected_genome_file + ".nin" not in all_files:
            genome_file_database = False

    # move on to next genome if database present
    if not ref_genome and genome_file_database is True:
        
        message = "\nGenome : %s blast database already exists" % genome
        logging.info(message)
        return genome_file_database
    
    # build database on existing genome fasta file if present
    if not ref_genome and expected_genome_file in all_files:

        message = "\nNOT detected database for : %s but %s file is present" \
                             " -> building database...\n" % (genome, expected_genome_file)
        logging.info(message)
        # make database
        make_blast_database(database_path, genome)
        # make sure to assign back this variable to TRUE
        genome_file_database = True
        return genome_file_database
    
    # download the genome as a tar file if both database and tar are missing
    if not ref_genome and genome_file_database is False and zip_file not in all_files:

        message = "\nGenome : %s blast database does not exist" \
              " -> downloading and decompressing the genome (it might take a couple of minutes)...\n" % genome
        logging.info(message)

        all_genomes = genomes_preprocessing.get_names(wdir, ref_species)
        accession_nr = [names[2] for names in all_genomes if genome in names][0]
        # download genome
        download_genome(genome, accession_nr, database_path)
        # make database
        make_blast_database(database_path, genome)

        # check if genome was downloaded and unpacked successfully
        genome_found = False
        for root, dirs, files in os.walk(database_path, topdown=False):
            for file in files:
                if file == genome + ".fasta":
                    genome_found = True

        # exit pipeline if fasta genome absent
        if not genome_found:
            genome_file_database = False
            message = "...FATAL_ERROR... : Genome : %s failed to download or decompress" \
                      " -> exiting the pipeline now...\n" % genome
            logging.info(message)
            return genome_file_database

        # validate that the database was generated
        all_files = []
        for root, dirs, files in os.walk(database_path, topdown=False):
            for f in files:
                all_files.append(f)

        for file in all_files:
            # look for indices
            if expected_genome_file + ".pal" not in all_files and expected_genome_file + ".nin" not in all_files:
                message = "\n...FATAL_ERROR... : Genome : %s database failed to build" \
                      " -> exiting the pipeline now...\n" % genome
                logging.info(message)
                genome_file_database = False
                return genome_file_database

        genome_file_database = True

        # fasta file and databases are present for the genome
        if genome_found and genome_file_database is True:
            genome_file_database = True
            message = "\nGenome : %s was downloaded and decompressed successfully" % genome
            logging.info(message)
            return genome_file_database

    # unpack and decompress reference genome
    if ref_genome:
        
        if expected_genome_file in all_files:

            message = "\nReference genome : %s is present\n" % genome
            logging.info(message)
            return True
            
        elif expected_genome_file not in all_files:

            message = "\nReference genome : %s does not exist" \
                              " -> downloading and decompressing it now...\n" % genome
            logging.info(message)
            all_genomes = genomes_preprocessing.get_names(wdir, ref_species, ref_genome=True)
            accession_nr = all_genomes[2]
            # download genome
            download_genome(genome, accession_nr, ref_genome_path)
            return True

    else:
        message = "Else statement when making blast database"
        logging.info(message)
        print("Else statement when making blast database")
        
    return genome_file_database


def download_genome(genome, accession_nr, database_path):
    """Downloads and unzips genomes using NCBI Datasets CLI, makes one fasta file"""

    # download genome
    start_time = time.time()
    filepath_1 = database_path + genome + ".zip"
    # need to exclude genomic cds because its also ".fna" which confuses the concatenation
    cmd1 = [pyinstaller_compatibility.resource_path("datasets"),
            "download", "genome", "accession", accession_nr, "--exclude-genomic-cds", "--filename", filepath_1]
    subprocess.call(cmd1, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))  # mute output and cautions

    # unzip all chromosomes into single file
    filepath_2 = database_path + genome + ".fasta"
    cmd2 = [pyinstaller_compatibility.resource_path("unzip"), "-pq", filepath_1, "cat", "*/*.fna", filepath_2]
    with open(filepath_2, "w") as outfile:
        subprocess.call(cmd2, stdout=outfile, stderr=open(os.devnull, 'wb'))  # redirect to file and mute cautions

    stop_time = time.time()
    message = "         -> Done : in %s min" % ((stop_time - start_time) / 60)
    logging.info(message)

    os.remove(filepath_1)


def make_blast_database(database_path, genome):
    """ Makes a blast database based on genome fasta file."""

    genome_file = genome + ".fasta"
    make_database_nucl = [pyinstaller_compatibility.resource_path("makeblastdb"),
                          "-in", database_path + genome_file, "-dbtype", "nucl"]
    make_database_prot = [pyinstaller_compatibility.resource_path("makeblastdb"),
                          "-in", database_path + genome_file, "-dbtype", "prot"]

    message = "\n                  Building blast database for genome : %s ..." % genome
    logging.info(message)
    subprocess.call(make_database_nucl, stdout=open(os.devnull, 'wb'))  # mute log info
    subprocess.call(make_database_prot, stdout=open(os.devnull, 'wb'))  # mute log info

