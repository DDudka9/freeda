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

Takes a MSA output, finds exons, clones exons and makes a final cds for protein
of a given genome.

"""

from Bio import AlignIO
from freeda import fasta_reader
from freeda import exon_finder
from freeda import cds_cloner
from freeda import input_extractor
from freeda import genomes_preprocessing
import operator
import logging
import re
import shutil
import glob

logging.basicConfig(level=logging.INFO, format="%(message)s")


def analyze_MSA(wdir, ref_species, MSA_path, gene, genome_name, ref_exons, expected_exons, all_genes_dict=None,
                final_excluded_species=None, subgroup=None):
    """Analyzes MSA per contig -> finds exons, clones them into cds"""

    # make a dictionary with exon number as key and sequences, names as values -> include microexons as empty lists
    microexons = input_extractor.check_microexons(wdir, gene, ref_species)

    final_exon_number = len(ref_exons)
    cloned_exons_overhangs = []
    # define pattern of all aligned fasta files
    pattern = "aligned*.fasta"
    # make a list of paths containing the aligned fasta files
    paths = glob.glob(MSA_path + "/" + pattern)

    for path in paths:
        # define contig name
        contig_name = re.search(r"(?<=aligned_).*$", path).group().replace("rev_comp_", "")
        contig_name = contig_name.rstrip(".fasta")
        message = "\n\n-------- %s ---------\n" % (contig_name)
        print(message)
        logging.info(message)
        # collect all sequences in that path
        seqs = collect_sequences(path)
        # define cds, locus for this contig, and gene
        cds_seq, locus_seq, gene_seq = index_positions(seqs)
        # find all exons in contig locus if no retrotransposition was detected
        exons, possible_retrotransposition, synteny, RETRO_score, synteny_score \
            = exon_finder.find_exons(ref_species, gene, cds_seq, locus_seq, gene_seq, contig_name, ref_exons,
                                     expected_exons, all_genes_dict)
        # skip this contig if possible retrotransposition event was detected
        # likelihood of false positive RETRO is more than 1 per 3 intronic exons
        # skip also contigs that are likely duplications or tandem repetitions
        if exons is None or (possible_retrotransposition is True and RETRO_score >= 0.4) or synteny_score < 0.5:
            continue
        # clone all exons WITH OVERHANGS
        cloned_exons_overhangs.append((contig_name, clone_exons_overhangs(seqs, exons)))
        # select only contigs carrying intronic exons
    preselected_exons_overhangs = preselect_exons_overhangs(cloned_exons_overhangs, expected_exons, microexons)
    # find contigs containing most intronic exons
    most_intronic_contigs = find_contigs_with_most_intronic_exons(preselected_exons_overhangs)
    
    # clone cds based on the most intronic contigs
    cloned_cds = cds_cloner.clone_cds(wdir, ref_species, preselected_exons_overhangs, most_intronic_contigs,
                                      gene, genome_name, final_exon_number, ref_exons, MSA_path)

    # check if final CDS is in frame (clone anyway)
    if (len(cloned_cds)-cloned_cds.count("-")) % 3 != 0:
        message = "\n !!!!! CDS is NOT in frame! !!!!!"
        print(message)
        logging.info(message)
        
    # remove dashes
    final_cds = cloned_cds.replace("-", "")  # THIS PIECE OF CODE IS OVERWRITTEN LATER
    
    # remove Ns
    if cloned_cds.count("N") > 0:
        side_note = "\n      N bases present in the cloned cds. They will be removed."
        print(side_note)
        logging.info(side_note)
    final_cds = cloned_cds.replace("N", "")
    
    # write a given gene cds into a file
    species_name = find_species_abbreviation(wdir, ref_species, gene, genome_name, final_cds, MSA_path,
                                             final_excluded_species, subgroup)
    file_cloned_cds(final_cds, gene, species_name)


def collect_sequences(path):
    """Collects sequences for a given path file (for each contig) in a list"""

    # use biopython to read alignment file
    alignment = AlignIO.read(path, "fasta")
    # collect sequences in the alignment
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]

    return seqs


def index_positions(seqs):
    """Gets cds, locus and gene sequences for this contig - each will be a dictionary with positions as values"""
    
    # for initial exon identification
    if len(seqs) == 4:
        a = seqs[2][1]
        b = seqs[3][1]

    # for later precise exon cloning
    else:
        a = seqs[1][1]
        b = seqs[2][1]
    
    cds_seq = {}
    locus_seq = {}
    gene_seq = {}

    nucleotides_read = 0 
    while nucleotides_read < len(seqs[0][1]):
        for position, nucleotide in enumerate(seqs[0][1]):
            cds_seq[position] = nucleotide
            nucleotides_read += 1

    # reset nucleotides read variable to 0
    nucleotides_read = 0
    while nucleotides_read < len(a):
        for position, nucleotide in enumerate(a):
            locus_seq[position] = nucleotide
            nucleotides_read += 1

    # reset nucleotides read variable to 0
    nucleotides_read = 0
    while nucleotides_read < len(b):
        for position, nucleotide in enumerate(b):
            gene_seq[position] = nucleotide
            nucleotides_read += 1

    return cds_seq, locus_seq, gene_seq


def clone_exons_overhangs(seqs, exons):
    """Gets the sequence of exons with overhangs"""

    cloned_exons_overhangs = []
    prefix = -50
    suffix = 50

    # read the exons dictionary to retrieve boundary positions
    for exon, features in exons.items():

        # unpack boundary tuple
        start, end, intron, exon_number, big_insertion = features

        # define locus sequence
        locus = seqs[2][1]
        gene = seqs[3][1]

        # clone exons
        locus_exon = locus[start+prefix:end+suffix]
        gene_exon = gene[start+prefix:end+suffix]
        cloned_exons_overhangs.append((locus_exon, gene_exon, intron, exon_number, big_insertion))

    return cloned_exons_overhangs


def preselect_exons_overhangs(cloned_exons_overhangs, expected_exons, microexons):
    """Preselects exons with overhangs"""

    intronic_exons = []
    preselected_exons_overhangs = {}
    sorted_exons = []
    for contig_exons in cloned_exons_overhangs:
        # unpack exons and their features per contig
        contig_name, exons = contig_exons
        
        # count each time intron changes in a locus
        missing_exon_detection = 0
        intr = False
        for exon, genomic, intron, exon_number, big_insertion in exons:
            if intron == intr:
                continue
            else:
                intr = intron
                missing_exon_detection += 1
                
        # log that possibly intron check was too strict for an exon in this contig
        if missing_exon_detection > 2:
            message = "\n       Possible intron missed by the code in contig: %s" % contig_name
            print(message)
            logging.info(message)

        intronic_exons.append([(final_exon_number, contig_name, exon, genomic) for
                               exon, genomic, intron, final_exon_number, big_insertion in exons
                               if intron is True and big_insertion is False])
            
        sorted_exons = sorted([entry for entry in intronic_exons if entry != []])

    all_exons = sorted(microexons + list(expected_exons))
    for number in range(1, len(all_exons) + 1):
        preselected_exons_overhangs[number] = []
        for exons in sorted_exons:
            for exon_nr, contig_name, exon, genomic in exons:
                if exon_nr == number:
                    preselected_exons_overhangs[number].append((contig_name, exon, genomic))
                    
    return preselected_exons_overhangs


def find_contigs_with_most_intronic_exons(preselected_exons_overhangs):
    """Finds contigs with most syntenic exons"""

    intronic_contigs = {}
    for exon, contigs in preselected_exons_overhangs.items():
        if contigs != []:
            for contig in contigs:
                try:
                    intronic_contigs[contig[0]] += 1
                except KeyError:
                    intronic_contigs[contig[0]] = 0
                    intronic_contigs[contig[0]] += 1

    most_intronic_contigs = sorted(intronic_contigs.items(), key=operator.itemgetter(1), reverse=True)
        
    return most_intronic_contigs


def find_species_abbreviation(wdir, ref_species, gene, genome_name, cloned_cds, MSA_path, final_excluded_species=None,
                              subgroup=None):
    """Finds species abbreviation."""

    # get all genome names and genomes dict with species abbreviations as keys
    all_genomes = genomes_preprocessing.get_names(wdir, ref_species, final_excluded_species, subgroup)

    header = [name[0] for name in all_genomes if name[1] == genome_name][0]
    filepath = MSA_path + "Cloned_cds"
    filename = header + "_" + gene + ".fasta"

    with open(filename, "w") as f:
        f.write(">" + header + "_" + gene + "_cds")
        f.write("\n" + cloned_cds)
    shutil.move(filename, filepath)

    return header


def file_cloned_cds(cloned_cds, gene, species_name):
    """Appends cloned cds to the end of the file with all cds."""

    file_name = gene + ".fasta"
    name = gene + "_" + species_name
    with open(file_name, "a+") as f:
        f.write(">" + name + "\n")
        f.write(cloned_cds + "\n")
    
    return file_name

