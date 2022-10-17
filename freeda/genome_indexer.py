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

