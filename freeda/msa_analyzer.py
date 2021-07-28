#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 18:13:46 2021

@author: damian

Takes a MAFFT output, finds exons, clones exons and makes a final cds for protein
of a given genome.

"""

from Bio import AlignIO
from freeda import fasta_reader
from freeda import exon_finder
from freeda import cds_cloner
import operator
import logging
import re
import shutil
import glob


def analyse_MSA(wdir, MSA_path, protein_name, genome_name, result_path, \
                Mm_exons, microexons, expected_exons):
    
    final_exon_number = len(Mm_exons)
    cloned_exons_overhangs = []
    # define pattern of all aligned fasta files
    pattern = "aligned*.fasta"
    # make a list of paths containins the aligned fasta files
    paths = glob.glob(MSA_path + "/" + pattern)
    for path in paths:
        # define contig name
        contig_name = re.search(r"(?<=aligned_).*$", path).group().replace("rev_comp_", "")
        contig_name = contig_name.rstrip(".fasta")
        message = "\n-------- %s ---------\n" % (contig_name)
        print(message)
        logging.info(message)
        # collect all sequences in that path
        seqs = collect_sequences(path)
        # define cds, locus for this contig, and gene
        cds, locus, gene = index_positions(seqs)
        # find all exons in contig locus if no retrotransposition was detected
        exons, possible_retrotransposition, synteny, RETRO_score, duplication_score \
            = exon_finder.find_exons(cds, locus, gene, seqs, contig_name, Mm_exons, expected_exons, microexons)
        # skip this contig if possible retrotransposition event was detected
        # likelihood of false positive RETRO is more than 1 per 3 intronic exons
        # skip also contigs that are likely duplications
        if (possible_retrotransposition == True and RETRO_score >= 0.4) or duplication_score < 0.5:
            continue
        # clone all exons WITH OVERHANGS
        cloned_exons_overhangs.append((contig_name, clone_exons_overhangs(seqs, exons)))
        # select only contigs carrying intronic exons
    preselected_exons_overhangs = preselect_exons_overhangs(cloned_exons_overhangs, final_exon_number, expected_exons)
    # find contigs containing most intronic exons
    most_intronic_contigs = find_contigs_with_most_intronic_exons(preselected_exons_overhangs)
    
    # clone cds based on the most intronic contigs
    cloned_cds = cds_cloner.clone_cds(preselected_exons_overhangs, most_intronic_contigs, \
                 protein_name, genome_name, final_exon_number, Mm_exons, MSA_path, microexons)

    # check if final CDS is in frame (clone anyway)
    if (len(cloned_cds)-cloned_cds.count("-")) % 3 != 0:
        message = "\n !!!!! CDS is NOT in frame! !!!!!"
        print(message)
        logging.info(message)
        
    # remove dashes
    final_cds = cloned_cds.replace("-", "")
    
    # remove Ns
    if cloned_cds.count("N") > 0:
        side_note = "\n      N bases present in the cloned cds. They will be removed."
        print(side_note)
        logging.info(side_note)
    final_cds = cloned_cds.replace("N", "")
    
    
    # write a given protein cds into a file
    species_name = write_cds(wdir, protein_name, genome_name, final_cds, MSA_path)
    file_cloned_cds(final_cds, protein_name, species_name)


def collect_sequences(path):
    # use biopython to read alignment file
    alignment = AlignIO.read(path, "fasta")
    # collect sequences in the alignment
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]
    return seqs


def index_positions(seqs):
    
    # for initial exon identification
    if len(seqs) == 4:
        a = seqs[2][1]
        b = seqs[3][1]
    # for later precise exon cloning
    else:
        a = seqs[1][1]
        b = seqs[2][1]
    
    cds = {}
    locus = {}
    gene = {}
    nucleotides_read = 0 
    while nucleotides_read < len(seqs[0][1]):
        for position, nucleotide in enumerate(seqs[0][1]):
            cds[position] = nucleotide
            nucleotides_read += 1
    # reset nucleotides read variable to 0
    nucleotides_read = 0
    while nucleotides_read < len(a):
        for position, nucleotide in enumerate(a):
            locus[position] = nucleotide
            nucleotides_read += 1
    # reset nucleotides read variable to 0
    nucleotides_read = 0
    while nucleotides_read < len(b):
        for position, nucleotide in enumerate(b):
            gene[position] = nucleotide
            nucleotides_read += 1
    return cds, locus, gene


def clone_exons_overhangs(seqs, exons): # works well
    cloned_exons_overhangs = []
    prefix = -50
    suffix = 50
    # read the exons dictionary to retrieve boundery positions 
    for exon, features in exons.items():
        # unpack boundery tuple
        start, end, introny, exon_number, big_insertion = features
        # define locus sequence
        locus = seqs[2][1]
        gene = seqs[3][1]
        # clone exons
        cloned_exons_overhangs.append((locus[start+prefix:end+suffix], \
                                       gene[start+prefix:end+suffix], \
                                        introny, exon_number, big_insertion))
    return cloned_exons_overhangs


def preselect_exons_overhangs(cloned_exons_overhangs, final_exon_number, expected_exons):
    intronic_exons = []
    preselected_exons_overhangs = {}
    sorted_exons = []
    for contig_exons in cloned_exons_overhangs:
        # unpack exons and their features per contig
        contig_name, exons = contig_exons
        
        # count each time introny changes in a locus
        missing_exon_detection = 0
        intr = False
        for exon, genomic, introny, exon_number, big_insertion in exons:
            if introny == intr:
                continue
            else:
                intr = introny
                missing_exon_detection += 1
                
        # log that possibly introny check was too strict for an exon in this contig
        if missing_exon_detection > 2:
            message = "\nPossible introny missed by the code in contig: %s" % contig_name
            print(message)
            logging.info(message)

        intronic_exons.append([(final_exon_number, contig_name, exon, genomic) for \
            exon, genomic, introny, final_exon_number, big_insertion in exons \
                             if introny == True and big_insertion == False]) 
            
        sorted_exons = sorted([entry for entry in intronic_exons if entry != []])
    
    # make a dictionary with exon number as key and sequences, names as values
    for number in expected_exons:

        preselected_exons_overhangs[number] = []
        for exons in sorted_exons:
            for exon_nr, contig_name, exon, genomic in exons:
                if exon_nr == number:
                    preselected_exons_overhangs[number].append((contig_name, exon, genomic))
                    
    return preselected_exons_overhangs


def find_contigs_with_most_intronic_exons(preselected_exons_overhangs): # works well
    intronic_contigs = {}
    for exon, contigs in preselected_exons_overhangs.items():
        if contigs != []:
            for contig in contigs:
                try:
                    intronic_contigs[contig[0]] += 1
                except KeyError:
                    intronic_contigs[contig[0]] = 0
                    intronic_contigs[contig[0]] += 1

    most_intronic_contigs = sorted(intronic_contigs.items(), key=operator.itemgetter(1), \
                                                                 reverse=True)
        
    return most_intronic_contigs


def write_cds(wdir, protein_name, genome_name, cloned_cds, MSA_path):
    filepath = MSA_path + "Cloned_cds"
    species_path = wdir + "species.txt"
    genomes_path = wdir + "genomes.txt"
    with open(species_path, "r") as s:
        species_list = s.readlines()
    with open(genomes_path, "r") as g:
        genomes_list = g.read().splitlines()
    
    species_index = genomes_list.index(genome_name + ".fasta")
    
    header = species_list[species_index].rstrip("\n")
    filename = header + "_" + protein_name + ".fasta"
    with open(filename, "w") as f:
        f.write(">" + header + "_" + protein_name + "_cds")
        f.write("\n" + cloned_cds)
    shutil.move(filename, filepath)
    return header


def file_cloned_cds(cloned_cds, protein_name, species_name):
    file_name = protein_name + ".fasta"
    name = protein_name + "_" + species_name
    with open(file_name, "a+") as f:
        f.write(">" + name + "\n")
        f.write(cloned_cds + "\n")
    
    return file_name