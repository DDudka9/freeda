#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 22:47:43 2021

@author: damian
"""

# REQUIRES INSTALLATION OF THE REFERENCE GENOME:
    # pyensembl install --species mouse --release 93
    # pyensembl install --species human --release 75

import pyensembl


def extract_input(wdir, original_species, protein):
    
    ensembl = pyensembl.EnsemblRelease(release=93, species = "mus_musculus")
    
    # find transcript
    transcript = extract_transcript(ensembl, protein)
    # find coding sequence
    coding_sequence = transcript.coding_sequence
    # find protein sequence
    protein_sequence = extract_protein(protein, transcript)
    # find gene sequence
    
    
    gene = pyensembl.Gene(gene_id, gene_name, contig, start, end, strand, biotype, genome)

    
    


def extract_protein(protein, transcript):
    """Extract protein sequnce based on transcript object"""
    
    # find protein id
    protein_id = transcript.protein_id
    # find proteins sequence
    protein_sequence = transcript.protein_sequence
    
    return protein_sequence
    
def extract:exons(transcript, longest_transcript_id):
    
    # create list od Exon objects
    exons = transcript.exons
    # get start and stop for each exon
    exon_intervals = transcript.exon_intervals
    # get sequence of each exon

    for exon in exon_intervals:
        print(exon[0])
        ensembl[exon[0]]
    

def extract_gene(ensembl, contig, transcript, gene_id):
    
    # list all genes in contig
    genes = ensembl.genes(contig)
    # create a Gene object
    gene = [gene for gene in genes if gene.gene_id == gene_id[0]][0]
    # get gene name
    gene_name = gene.gene_name

    
def extract_transcript(ensembl, protein):
    """Extracts transcript for a given protein by creating transcript object"""
    
    all_transcripts_ids = ensembl.transcript_ids_of_gene_name(protein)
    # pick the longest transcript
    transcript_id = sorted(all_transcripts)[0]
    # define biotype (FREEDA deals only with protein coding sequences)
    biotype = "Protein coding"
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
    transcript = pyensembl.Transcript(transcript_id, transcript_name, contig, start, end, strand, biotype,
                               gene_id, ensembl, support_level=None)
    # check if given transcript has an annotated start codon
    start_codon_present = transcript.contains_start_codon
    # check if given transcript has an annotated stop codon
    stop_codon_present = transcript.contains_stop_codon
    
    if start_codon_present and stop_codon_present:
        return transcript
    else:
        print("Chosen transcript for protein: %s does not have START and STOP annotated" % protein)
        
        
        
        
        
        
        
        
        
    
      
    
    