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

def run_blast(wdir, original_species):
    genomes_file_dir = wdir + "genomes.txt" 
    proteins_file_dir = wdir + "proteins.txt" 
    
    database_path = wdir + "Genomes/"
    query_path = wdir + "Blast_input/"
    output_path = wdir + "Blast_output/"
    form = "6 qseqid means sseqid means qstart means qend means sstart means send means evalue means bitscore length means pident means mismatch means gapopen means qlen means slen means"
    
    genomes = [genome.rstrip("\n") for genome in open(genomes_file_dir, "r").readlines()]
    proteins = open(proteins_file_dir, "r").readlines()
    
    # record all the existing files in the blast output directory
    all_files = set(os.listdir(output_path))
    
    for genome in genomes:
        for protein in proteins:
            pattern = protein + "_" + genome
            
            # to avoid duplications
            if [string for string in all_files if pattern in string] == []:
                
                protein = protein.lstrip("\n").rstrip("\n")
                database = database_path + genome
                query = query_path + protein + "_" + original_species + "_uniprot.fasta"
                output = output_path + protein + "_" + genome + ".txt"
                to_blast = ["tblastn", "-db", database, "-query", query, "-out", output, "-outfmt", form, "-num_threads", "5"]
                subprocess.call(to_blast)
            
            else:
                print("Protein: %s from genome: %s -> was already blasted." % (protein, genome))
    
    print("tblastn txt files have been generated.")
    
    return output_path