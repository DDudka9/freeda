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
#import glob
import shutil
wdir = os.getcwd() + "/"
original_species = "Mm"
protein = "Cenpt"

rules = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
         "Y": "R", "R": "Y", "W": "W", "S": "S", "K": "M", "M": "K",
         "D": "H", "H": "D", "V": "B", "B": "V", "X": "X"}

def extract_input(wdir, original_species, genome, protein):
    """Extracts all input sequences required by FREEDA from the indicated reference genome (original_species)"""
    
    blast_input_path = wdir + "Blast_input/"
    coding_sequence_input_path = wdir + "Coding_sequences/"
    gene_input_path = wdir + "Genes/"
    exons_input_path = wdir + "Exons/"
    
    # default is mouse
    mouse_names = {"Mm", "mouse", "Mus musculus", "mus musculus"}
    human_names = {"Hs", "human", "Homo sapiens", "homo sapiens"}
    if original_species in mouse_names:
        #header = "musculus"
        species = "mus musculus"
        release = 93
        reference_genome_contigs_file = "MUSCULUS_genome_contigs.txt"
        reference_genome_name = "MUSCULUS_genome"
    
    elif original_species in human_names:
        #header = "sapiens"
        species = "homo sapiens"
        release = 75
        reference_genome_contigs_file = "SAPIENS_genome_contigs.txt"
        reference_genome_name = "SAPIENS_genome"
    
    # make sure original species genome (reference genome) is present
    reference_genomes_path = wdir + "Reference_genomes/"
    # find reference genome directory based on the header variable
    #reference_genome_dir = [filename for filename in glob.glob(reference_genomes_path + "*") if header in filename.lower()]
    # extract name of the genome
    #reference_genome_name = reference_genome_dir[0].split("/")[-1].split(".")[0]
    # check if reference genome is present -> exit 
    reference_genome_present = check_genome_present(reference_genomes_path, reference_genome_name, reference_genome=True)

    # get assembly database
    ensembl = pyensembl.EnsemblRelease(release, species)
    # define biotype (FREEDA deals only with protein coding sequences)
    biotype = "Protein coding"
    
    # find coding sequence
    transcript, longest_transcript_id, gene_id, contig, strand, UTR_5, UTR_3 = extract_cds(ensembl, 
                original_species, coding_sequence_input_path, protein, biotype)
    # find coding sequence
    #extract_cds(ensembl, original_species, coding_sequence_input_path, protein, biotype)
    # find protein sequence
    extract_protein(original_species, blast_input_path, protein, strand, transcript)
    # find gene sequence
    reference_contigs_dict =  extract_gene(original_species, gene_input_path, 
                                ensembl, contig, strand, gene_id, reference_genomes_path, 
                                reference_genome_name, reference_genome_contigs_file, protein)
    # find exons sequence
    extract_exons(protein, exons_input_path, reference_genomes_path, reference_genome_name, 
                  contig, strand, transcript, reference_contigs_dict, UTR_5, UTR_3)
    
    return reference_genome_present
    
    
def extract_protein(original_species, blast_input_path, protein, strand, transcript):
    """Extracts protein sequence based on Transcript object"""
    
    # find proteins sequence
    protein_sequence = transcript.protein_sequence
    # get and save the sequence
    parse_sequence(blast_input_path, protein_sequence, protein, strand = None, sequence_type = "protein")
    
    
def extract_exons(protein, exons_input_path, reference_genomes_path, reference_genome_name, 
                  contig, strand, transcript, reference_contigs_dict, UTR_5, UTR_3):
    
    # TEST "+" STRAND SEQUENCES !!!!
    # FIGURE OUT HOW TO DEAL WITH UTRs SPANNING MULTIPLE EXONS !!!
    
    
    """Extracts exoms sequence based on Transcript object"""
        
    # get start and stop for each exon
    exon_intervals = transcript.exon_intervals
   
    # get sequence of each exon
    exons_filename = protein + "_" + original_species + "_exons.fasta"
    exons_file = open(exons_input_path + exons_filename, "w")
    
    number = 0
    
    for exon in exon_intervals:
        number += 1
        start = exon[0]
        stop = exon[1]
        header = protein + "_" + original_species + "_exon" + str(number)
        exon_bed_filename = header + ".bed"
        exon_fasta_filename = header + ".fasta"
        with open(exon_bed_filename, "w") as b:
            b.write(reference_contigs_dict[contig] + "\t")
            
            if strand == "-":
                start = start - 1

            # TEST "+" STRAND SEQUENCES !!!!
            
            b.write(str(start) + "\t")
            b.write(str(stop))
        
        bed_object = pybedtools.BedTool(exon_bed_filename)
        bed_object = bed_object.sequence(fi = reference_genomes_path + \
                            reference_genome_name + ".fasta")
        exon_fasta_sequence = bed_object.save_seqs(exon_fasta_filename)
        
        with open(exon_fasta_filename, "r") as f:
            file = f.readlines()
            exon_fasta_sequence = file[1].rstrip("\n").upper()
               
            # write into file with all exons
            exons_file.write(">" + header + "\n")
            if strand == "-":
                complement = ""
                for base in exon_fasta_sequence:
                    complement = rules[base] + complement
                    exon_fasta_sequence = complement
            
            if number == 1:
                # already in correct orientation
                exon_fasta_sequence_trimmed = exon_fasta_sequence.replace(UTR_5, "")
                exon_fasta_sequence = exon_fasta_sequence_trimmed
                
            if number == len(exon_intervals):
                # already in correct orientation
                exon_fasta_sequence_trimmed = exon_fasta_sequence.replace(UTR_3, "")
                exon_fasta_sequence = exon_fasta_sequence_trimmed
                
            exons_file.write(exon_fasta_sequence + "\n")
        
        os.remove(exon_bed_filename)
        os.remove(exon_fasta_filename)

    exons_file.close()
    

def extract_gene(original_species, gene_input_path, ensembl, contig, strand, gene_id, 
        reference_genomes_path, reference_genome_name, reference_genome_contigs_file, protein):
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
    gene_fasta_sequence = bed_object.sequence(fi = reference_genomes_path + \
                            reference_genome_name + ".fasta")
        
    os.remove(gene_bed_filename)
        
    parse_sequence(gene_input_path, gene_fasta_sequence, protein, strand, sequence_type="gene")
    
    return reference_contigs_dict

    #print(open(gene_fasta_file.seqfn).read())
    #gene_sequence = gene_fasta_sequence_file.save_seqs(gene_fasta_filename)
    

def parse_sequence(output_path, fasta_sequence, protein, strand, sequence_type):
    """Parses either a pybedtools.bedtool.BedTool object or string into a file"""
    
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
    
    # TEST LONGEST TRANSCRIPT FINDING ON PROTEINS WITH DIFFERENT TRANSCRIPT LENGTHS
    
    all_transcripts_ids = ensembl.transcript_ids_of_gene_name(protein)
    all_transcripts_dict = {}
    
    # find longest transcript
    for t in all_transcripts_ids:
        
        all_transcripts_dict[t] = []
    
        # find gene id
        gene_id = ensembl.gene_ids_of_gene_name(protein)
        all_transcripts_dict[t].append(gene_id)
        
        # find contig
        contig = ensembl.locus_of_transcript_id(t).contig
        all_transcripts_dict[t].append(contig)
        
        # find strand
        strand = ensembl.locus_of_transcript_id(t).strand
        all_transcripts_dict[t].append(strand)
        
        # find contig name
        transcript_name = ensembl.transcript_name_of_transcript_id(t)
        all_transcripts_dict[t].append(transcript_name)
        
        # find start of the transcript (not START codon)
        start = ensembl.locus_of_transcript_id(t).start
        all_transcripts_dict[t].append(start)
        
        # find end of the transcript (not STOP codon)
        end = ensembl.locus_of_transcript_id(t).end
        all_transcripts_dict[t].append(end)
        
        # find transcript
        transcript = pyensembl.Transcript(t, transcript_name, 
            contig, start, end, strand, biotype, gene_id, ensembl, support_level=None)
        
        # try to get cds -> if "start_codon" does not exist -> KeyError and then ValueError occur
        try:
            cds_sequence = transcript.coding_sequence
        except ValueError:
            cds_sequence = ""
            all_transcripts_dict[t].append(cds_sequence)
        finally:
            all_transcripts_dict[t].append(cds_sequence)
        
        # get length of cds
        length = len(cds_sequence)
        all_transcripts_dict[t].append(length)
        
    # find transcript with the longest coding sequence
    longest_transcript_id = None
    length = 0
    for t, features in all_transcripts_dict.items():
        if features[-1] > length:
            longest_transcript_id = t
            length = features[-1]
    
    # unpack features of the longest transcript
    gene_id, contig, strand, transcript_name, start, end, cds_sequence, length = all_transcripts_dict[longest_transcript_id]
    # rebuild Transcript object based on the longest transcript
    transcript = pyensembl.Transcript(longest_transcript_id, transcript_name, 
            contig, start, end, strand, biotype, gene_id, ensembl, support_level=None)
    
    # check if given transcript has an annotated start codon
    start_codon_present = transcript.contains_start_codon
    # get 5'UTR sequence
    UTR_5 = transcript.five_prime_utr_sequence.upper()
    # check if given transcript has an annotated stop codon
    stop_codon_present = transcript.contains_stop_codon
    # get 3'UTR sequence
    UTR_3 = transcript.three_prime_utr_sequence.upper()
    # get and save the sequence
    parse_sequence(coding_sequence_input_path, cds_sequence, protein, 
                                              strand, sequence_type="cds")
    
    if start_codon_present and stop_codon_present:
        return transcript, longest_transcript_id, gene_id, contig, strand, UTR_5, UTR_3
    else:
        print("Chosen transcript for protein: %s does not have START and STOP annotated" % protein)
        
        
        
        
        
        
        
        
        
    
      
    
    