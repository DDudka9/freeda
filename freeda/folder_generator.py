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

Generates the initial folders before input can be extracted

"""

import os


def generate_basic_folders(wdir):
    """Checks if folders for input are present in working directory, generates if not"""

    path_to_ref_genomes = wdir + "Reference_genomes/"
    path_to_genomes = wdir + "Genomes/"
    path_to_blast_input = wdir + "Blast_input/"
    path_to_blast_output = wdir + "Blast_output/"
    path_to_cds = wdir + "Coding_sequences/"
    path_to_genes = wdir + "Genes/"
    path_to_exons = wdir + "Exons/"
    path_to_structures = wdir + "Structures/"

    if not os.path.isdir(path_to_ref_genomes):
        os.makedirs(path_to_ref_genomes)
    if not os.path.isdir(path_to_genomes):
        os.makedirs(path_to_genomes)
    if not os.path.isdir(path_to_blast_input):
        os.makedirs(path_to_blast_input)
    if not os.path.isdir(path_to_blast_output):
        os.makedirs(path_to_blast_output)
    if not os.path.isdir(path_to_cds):
        os.makedirs(path_to_cds)
    if not os.path.isdir(path_to_exons):
        os.makedirs(path_to_exons)
    if not os.path.isdir(path_to_genes):
        os.makedirs(path_to_genes)
    if not os.path.isdir(path_to_structures):
        os.makedirs(path_to_structures)


def generate_contig_folders(result_path, all_genes, all_genomes):
    """Generates folders for contigs"""
    
    # make gene folders
    for gene in all_genes:
        os.makedirs(result_path + str(gene))

    # name genome folder for each gene folder
    all_gene_folders = [folder for folder in os.listdir(result_path) if not folder.startswith(".")]

    for gene_folder in all_gene_folders:
        for genome in all_genomes:
            os.makedirs(result_path + str(gene_folder) + "/" + genome)
            # make Bed folder for each genome folder per gene
            os.makedirs(result_path + str(gene_folder) + "/" + genome + "/" + "Bed/above_threshold")
            # make Fasta folder for each genome folder per gene
            os.makedirs(result_path + str(gene_folder) + "/" + genome + "/" + "Fasta/above_threshold")
            # make MSA folder for each genome folder per gene
            os.makedirs(result_path + str(gene_folder) + "/" + genome + "/" + "Fasta/above_threshold/MSA")
            # make Cloned_exons folder
            os.makedirs(result_path + str(gene_folder) + "/" + genome + "/" + "Fasta/above_threshold/MSA/Cloned_cds")
            # make Duplicated_exons folder
            os.makedirs(result_path + str(gene_folder) + "/" + genome + "/" + "Fasta/above_threshold/MSA/Duplicated_exons")
            # make Single_exon_MSA folder
            os.makedirs(result_path + str(gene_folder) + "/" + genome + "/" + "Fasta/above_threshold/MSA/Single_exons_MSA")

    # make folders for output files
    result_path = result_path.replace("Raw_data/", "Results/")
    os.makedirs(result_path + "Nucleotide_alignments")
    os.makedirs(result_path + "Protein_alignments")
    os.makedirs(result_path + "Graphs")
    os.makedirs(result_path + "Results_sheet")
    os.makedirs(result_path + "Gene_trees")
    os.makedirs(result_path + "Structures")
    os.makedirs(result_path + "Blast_output")

    print("\nFolders have been generated.\n")