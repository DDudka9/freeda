#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 22:47:43 2021

@author: damian
"""

# REQUIRES INSTALLATION OF THE REFERENCE GENOME:
    # pyensembl install --species mouse --release 93
    # pyensembl install --species human --release 75

from freeda.tblastn import check_genome_present
import pyensembl
import pybedtools

import os
import glob
wdir = os.getcwd() + "/"
original_species = "Mm"
protein = "Cenpt"

def extract_input(wdir, original_species, genome, protein):
    """Extracts all input sequences required by FREEDA from the indicated reference genome (original_species)"""
    
    blast_input_path = wdir + "Blast_input/"
    coding_sequence_input_path = wdir + "Coding_sequences/"
    gene_input_path = wdir + "Genes/"
    exons_input_path = wdir + "Exons/"
    
    if original_species == "Hs":
        header = "sapiens"
        species = "homo sapiens"
        release = 75
        reference_genome_contigs_file = "SAPIENS_genome_contigs.txt"
    
    # default is mouse
    else:
        header = "musculus"
        species = "mus musculus"
        release = 93
        reference_genome_contigs_file = "MUSCULUS_genome_contigs.txt"
    
    # make sure original species genome (reference genome) is present
    reference_genomes_path = wdir + "Reference_genomes/"
    # find reference genome directory based on the header variable
    reference_genome_dir = [filename for filename in glob.glob(reference_genomes_path + "*") if header in filename.lower()]
    # extract name of the genome
    reference_genome_name = reference_genome_dir[0].split("/")[-1].split(".")[0]
    # check if reference genome is present -> exit 
    reference_genome_present = check_genome_present(reference_genomes_path, reference_genome_name, reference_genome=True)

    # get assembly database
    ensembl = pyensembl.EnsemblRelease(release, species)
    # define biotype (FREEDA deals only with protein coding sequences)
    biotype = "Protein coding"
    
    # find coding sequence
    transcript, gene_id, contig, strand = extract_cds(ensembl, original_species, 
                                    coding_sequence_input_path, protein, biotype)
    # find protein sequence
    extract_protein(original_species, blast_input_path, protein, strand, transcript)
    # find gene sequence
    extract_gene(original_species, gene_input_path, ensembl, contig, strand, gene_id, 
                 reference_genomes_path, reference_genome_name, protein)
    
    
def extract_protein(original_species, blast_input_path, protein, strand, transcript):
    """Extracts protein sequence based on Transcript object"""
    
    # find proteins sequence
    protein_sequence = transcript.protein_sequence
    # get and save the sequence
    parse_sequence(blast_input_path, protein_sequence, protein, strand = None, sequence_type = "protein")
    
    
def extract:exons(transcript, longest_transcript_id):
    """Extracts exoms sequence based on Transcript object"""
    
    # create list od Exon objects
    exons = transcript.exons
    # get start and stop for each exon
    exon_intervals = transcript.exon_intervals
    # get sequence of each exon

    for exon in exon_intervals:
        print(exon[0])
        ensembl[exon[0]]
    

def extract_gene(original_species, gene_input_path, ensembl, contig, strand, gene_id, 
                 reference_genomes_path, reference_genome_name, protein):
    """Extracts gene sequence based on Genome object"""
    
    # list all genes in contig
    genes = ensembl.genes(contig)
    # create a Gene object
    gene = [gene for gene in genes if gene.gene_id == gene_id[0]][0]
    # get gene name
    gene_name = gene.gene_name
    # get gene starting position
    start = gene.start
    # get gene end position
    end = gene.end
    # get gene sequence
    
    # make a bed and fasta file for gene (add underscore to differenciate from other handles)
    gene_bed_filename = "_" + gene_name + ".bed"
    gene_fasta_filename = "_" + gene_name + ".fasta"
    
    # match contigs named between pyensembl Genome object and GeneBank assembly
    reference_contigs_dict = {}
    with open(reference_genome_contigs_file, "r") as f:
        file = f.readlines()
        for line in file:
            l = line.split("\t")
            reference_contigs_dict[l[0]] = l[4]
        
    matched_contig = reference_contigs_dict[str(contig)]
    
    # make a bedtool file
    with open(gene_bed_filename, "w") as b:
        b.write(str(matched_contig) + "\t" + str(start) + "\t" + str(end))
    
    # parse that bedtool file using BedTools object
    bed_object = pybedtools.BedTool(gene_bed_filename)
    gene_fasta_sequence_file = bed_object.sequence(fi = reference_genomes_path + \
                            reference_genome_name + ".fasta")
    os.remove(gene_bed_filename)
        
    parse_sequence(gene_input_path, gene_fasta_sequence_file, protein, strand, sequence_type="gene")
    
        
    #print(open(gene_fasta_file.seqfn).read())
    #gene_sequence = gene_fasta_sequence_file.save_seqs(gene_fasta_filename)
    

def parse_sequence(output_path, fasta_sequence, protein, strand, sequence_type):
    """Parses either a pybedtools.bedtool.BedTool object or string into a file"""
    
    rules = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
             "Y": "R", "R": "Y", "W": "W", "S": "S", "K": "M", "M": "K",
             "D": "H", "H": "D", "V": "B", "B": "V", "X": "X"}
    
    if isinstance(fasta_sequence, pybedtools.bedtool.BedTool):
        file = open(fasta_sequence.seqfn).readlines()
        header = ">" + protein + "_" + original_species + "_" + sequence_type
        sequence = file[1].upper().rstrip("\n")
    
    if not isinstance(fasta_sequence, pybedtools.bedtool.BedTool):
        header = ">" + protein + "_" + original_species + "_" + sequence_type
        sequence = fasta_sequence
        
    if sequence_type == "gene" and strand == "-":
        complement = ""
        for base in sequence:
            complement = rules[base] + complement
        sequence = complement
    
    filename = protein + "_" + original_species + "_" + sequence_type
    with open(output_path + filename + ".fasta", "w") as f:
        f.write(header + "\n")
        f.write(sequence)

    
def extract_cds(ensembl, original_species, coding_sequence_input_path, protein, biotype):
    """Extracts coding sequence by creating Transcript object"""
    
    all_transcripts_ids = ensembl.transcript_ids_of_gene_name(protein)
    # pick the longest transcript
    longest_transcript_id = sorted(all_transcripts_ids)[0]
    # find gene id
    gene_id = ensembl.gene_ids_of_gene_name(protein)
    # find contig
    contig = ensembl.locus_of_transcript_id(longest_transcript_id).contig
    # find strand
    strand = ensembl.locus_of_transcript_id(longest_transcript_id).strand
    # find contig name
    transcript_name = ensembl.transcript_name_of_transcript_id(longest_transcript_id)
    # find start
    start = ensembl.locus_of_transcript_id(longest_transcript_id).start
    # find end
    end = ensembl.locus_of_transcript_id(longest_transcript_id).end
    # find transcript
    transcript = pyensembl.Transcript(longest_transcript_id, transcript_name, contig, start, end, strand, biotype,
                               gene_id, ensembl, support_level=None)
    # check if given transcript has an annotated start codon
    start_codon_present = transcript.contains_start_codon
    # check if given transcript has an annotated stop codon
    stop_codon_present = transcript.contains_stop_codon
    
    # get cds
    cds_sequence = transcript.coding_sequence
    # get and save the sequence
    parse_sequence(coding_sequence_input_path, cds_sequence, protein, strand, sequence_type="cds")
    
    if start_codon_present and stop_codon_present:
        return transcript, gene_id, contig, strand
    else:
        print("Chosen transcript for protein: %s does not have START and STOP annotated" % protein)
        
        
        
        
        
        
        
        
        
    
      
    
    