#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:10:33 2021

@author: Damian Dudka - damiandudka0@gmail.com

Indexes a genome of interest using a BioPython module SeqIO

"""

from Bio import SeqIO
import logging
import glob


def index_genome_database(wdir, genome_name):
    """Indexes genome databases"""
    
    genomes_dir = wdir + "Genomes/"
    
    name = glob.glob(genomes_dir + genome_name + ".fasta")
    index = glob.glob(genomes_dir + genome_name + ".idx")

    if index == []:
        # it will make a new index
        message = "\nIndexing genome: " + genome_name + "\n"
        print(message)
        logging.info(message)
        genome_index = SeqIO.index_db(genomes_dir + genome_name + ".idx", name[0], "fasta")
        message = "\nGenome: " + genome_name + " has been indexed.\n"
        print(message)
        logging.info(message)
        
    else:
        # it will find the index and put it in the genome_index variable
        genome_index = SeqIO.index_db(genomes_dir + genome_name + ".idx", name[0], "fasta")
        message = "\nGenome: " + genome_name + " index already exists.\n"
        print(message)
        logging.info(message)
        
    return genome_index

