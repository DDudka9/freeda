#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:21:53 2021

@author: damian

Runs blast using NCBI tblastn. 
Requires protein database built using queried genomes:
https://www.ncbi.nlm.nih.gov/books/NBK279688/
Fasta genome names need to be stored in "genomes.txt" file one per line using one underscore:
    
    xxxxx_xxxxx.fasta
    yyyyy_yyyyy.fasta
    
Protein names need to be stored in "proteins.txt" file one per line:
    
    aaaaa
    bbbbb

Outputs tabulated text files per protein per genomes.

"""

import subprocess
import os
import glob
import tarfile
import shutil
import gzip


def run_blast(wdir, original_species):
    """Runs tblastn based on NCBI makedatabase routine."""

    genomes_file_dir = wdir + "genomes.txt"
    proteins_file_dir = wdir + "proteins.txt"

    database_path = wdir + "Genomes/"
    query_path = wdir + "Blast_input/"
    output_path = wdir + "Blast_output/"
    form = "6 qseqid means sseqid means qstart means qend means sstart means send means evalue "
    "means bitscore length means pident means mismatch means gapopen means qlen means slen means"

    genomes = [genome.rstrip("\n") for genome in open(genomes_file_dir, "r").readlines()]
    proteins = open(proteins_file_dir, "r").readlines()

    # clear Blast_output folder
    all_old_blast_output_files = glob.glob(os.path.join(output_path, "*.txt"))
    for file in all_old_blast_output_files:
        os.remove(file)

    # make sure database for each genome already exists or is successfully built
    for genome in genomes:
        genome_file_database = check_genome_present(database_path, genome.rstrip(".fasta"))
        if not genome_file_database:
            return None

    # perform blast
    for genome in genomes:

        for protein in proteins:
            protein = protein.lstrip("\n").rstrip("\n")
            database = database_path + genome
            query = query_path + protein + "_" + original_species + "_protein.fasta"
            output = output_path + protein + "_" + genome + ".txt"
            to_blast = ["tblastn", "-db", database, "-query", query, "-out", output, "-outfmt", form, "-num_threads",
                        "8"]
            print("\nPerforming tblast for protein: %s from genome: %s\n" % (protein, genome))
            subprocess.call(to_blast)

    print("\ntblastn txt files have been generated.")

    return output_path


def check_genome_present(database_path, genome, reference_genome=False):
    """Checks if a given genome is present. Unpacks and unzips genomes downloaded from NCBI Assembly.
    
    Non-ncbi assemblies must be prepared as ".fasta" files conform with "genomes.txt" names.
    
    It also looks for reference genome if key-only argument reference_genome is invoked."""

    genome_file_database = True

    if reference_genome:
        reference_genome_path = database_path

    # define expected files
    tar_file = genome + ".tar"
    expected_genome_file = genome + ".fasta"
    expected_database_extensions = {".00.phr",
                                    ".00.pin",
                                    ".00.psq",
                                    ".01.phr",
                                    ".01.pin",
                                    ".01.psq",
                                    ".02.phr",
                                    ".02.pin",
                                    ".02.psq",
                                    ".nhr",
                                    ".nin",
                                    ".nsq",
                                    ".pal"}

    # get info on all tar files available
    all_files = []
    for root, dirs, files in os.walk(database_path, topdown=False):
        for f in files:
            all_files.append(f)

    # check if genome database is present
    for file in all_files:
        for extension in expected_database_extensions:
            if expected_genome_file + extension not in all_files:
                genome_file_database = False

    # move on to next genome if database present
    if reference_genome is False and genome_file_database is True:

        print("\nGenome : %s blast database already exists" % genome)

        return genome_file_database

    # build datbase on existing genome fasta file if present
    elif reference_genome is False and expected_genome_file in all_files:

        print("\nNOT detected database for : %s but %s file is present"
              " -> building database...\n" % (genome, expected_genome_file))

        # make database
        make_blast_database(database_path, genome)
        # make sure to tick back this parameter to TRUE
        genome_file_database = True

        return genome_file_database

    # exit pipeline if both database and tar are missing
    elif reference_genome is False and genome_file_database is False and tar_file not in all_files:

        print("\nGenome : %s blast database does not exists and archive file is missing"
              " -> exiting the pipeline now...\n" % genome)

        return genome_file_database

    # unpack and unzip if tar file present
    elif reference_genome is False and genome_file_database is False and tar_file in all_files:

        print("\nNOT detected database for : %s -> building database...\n " % genome)

        # unpack and decompress the genome into a fasta file
        genome_to_unzip = database_path + "/" + genome + ".tar"
        unzip_genome(database_path, genome, genome_to_unzip)
        # remove the tar file
        os.remove(database_path + tar_file)
        # make database
        make_blast_database(database_path, genome)
        # make sure to tick back this parameter to TRUE
        genome_file_database = True

    # unpack and decompress reference genome
    elif reference_genome is True:

        if expected_genome_file in all_files:
            print("\nReference genome : %s is present\n" % genome)
            return True

        elif expected_genome_file not in all_files and tar_file in all_files:

            print("\nReference genome : %s does not exists but archive file is present"
                  " -> unpacking and decompressing...\n" % genome)

            # unpack and decompress the genome into a fasta file
            genome_to_unzip = reference_genome_path + "/" + genome + ".tar"
            unzip_genome(reference_genome_path, genome, genome_to_unzip)
            # remove the tar file
            os.remove(reference_genome_path + tar_file)
            return True

        elif expected_genome_file not in all_files and tar_file not in all_files:
            print("\nReference genome : %s does not exists and archive file is missing"
                  " -> exiting the pipeline now...\n" % genome)
            return False

    else:
        print("Else statement when making blast database")

    return genome_file_database


def unzip_genome(database_path, genome, genome_to_unzip):
    """Unzips tar archive files using gzip and shutil modules"""

    tar = tarfile.open(genome_to_unzip)
    for tarinfo in tar:
        if os.path.splitext(tarinfo.name)[1] == ".gz":
            source_directory = database_path + tarinfo.name
            filename = tarinfo.name.split("/")[1]

    tar.extractall(path=database_path, members=find_packed_genome(tar, genome))
    tar.close()

    shutil.move(source_directory, database_path)
    # remove the now empty directory
    moved_gz_filename = source_directory.split("/")[-1]
    empty_directory = source_directory.rstrip(moved_gz_filename)
    os.rmdir(empty_directory)
    # open the genome content in text mode ("rt")
    with gzip.open(database_path + filename, 'rt') as f_in:
        print("\n        Detected compressed genome file for : %s -> decompressing...\n " % genome)

        # write it into a fasta file
        with open(database_path + genome + '.fasta', 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)

    # remove compressed file
    os.remove(database_path + filename)


def find_packed_genome(members, genome):
    """Finds a packed genome in the archive file using tarfile module."""

    for tarinfo in members:
        if os.path.splitext(tarinfo.name)[1] == ".gz":
            print("\n    Detected archive genome file for : %s -> unpacking...\n " % genome)

            yield tarinfo


def make_blast_database(database_path, genome):
    """ Makes a blast database based on genome fasta file."""

    genome_file = genome + ".fasta"
    make_database_nucl = ["makeblastdb", "-in", database_path + genome_file, "-dbtype", "nucl"]
    subprocess.call(make_database_nucl)
    make_database_prot = ["makeblastdb", "-in", database_path + genome_file, "-dbtype", "prot"]
    subprocess.call(make_database_prot)


"""
            # to avoid duplications
            if [string for string in all_files if pattern in string] == []:
                
                protein = protein.lstrip("\n").rstrip("\n")
                database = database_path + genome
                query = query_path + protein + "_" + original_species + "_protein.fasta"
                output = output_path + protein + "_" + genome + ".txt"
                to_blast = ["tblastn", "-db", database, "-query", query, "-out", output, "-outfmt", form, "-num_threads", "5"]
                subprocess.call(to_blast)
            
            else:
                print("\nProtein: %s from genome: %s -> was already blasted." % (protein, genome))    
    
"""
