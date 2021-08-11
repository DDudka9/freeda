#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:04:31 2021

@author: damian
"""

import os


def generate_folders(result_path):
    # make protein folders
    for protein in open("proteins.txt", "r").readlines():
        os.makedirs(result_path + str(protein.rstrip("\n")))

    # name genome folder for each protein folder
    all_protein_folders = os.listdir(result_path)
    for protein_folder in all_protein_folders:
        for genome in open("genomes.txt", "r"):
            os.makedirs(result_path + str(protein_folder)
                        + "/" + str(genome.rstrip(".fasta\n")))
            # make Bed folder for each genome folder per protein
            os.makedirs(result_path + str(protein_folder)
                        + "/" + str(genome.rstrip(".fasta\n"))
                        + "/" + "Bed/above_threshold")
            # make Fasta folder for each genome folder per protein
            os.makedirs(result_path + str(protein_folder)
                        + "/" + str(genome.rstrip(".fasta\n"))
                        + "/" + "Fasta/above_threshold")
            # make MSA folder for each genome folder per protein
            os.makedirs(result_path + str(protein_folder)
                        + "/" + str(genome.rstrip(".fasta\n"))
                        + "/" + "Fasta/above_threshold/MSA")
            # make Cloned_exons folder
            os.makedirs(result_path + str(protein_folder)
                        + "/" + str(genome.rstrip(".fasta\n"))
                        + "/" + "Fasta/above_threshold/MSA/Cloned_cds")
            # make Duplicated_exons folder
            os.makedirs(result_path + str(protein_folder)
                        + "/" + str(genome.rstrip(".fasta\n"))
                        + "/" + "Fasta/above_threshold/MSA/Duplicated_exons")
            # make Single_exon_MSA folder
            os.makedirs(result_path + str(protein_folder)
                        + "/" + str(genome.rstrip(".fasta\n"))
                        + "/" + "Fasta/above_threshold/MSA/Single_exons_MSA")

    print("Folders have been generated.\n")
