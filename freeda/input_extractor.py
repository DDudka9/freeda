#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 22:47:43 2021

@author: damian
"""

# REQUIRES INSTALLATION OF THE REFERENCE GENOME:
    # pyensembl install --species mouse --release 100
    # pyensembl install --species human --release 100

from freeda.tblastn import check_genome_present
from freeda import reference_genome_dict_generator
import pyensembl
import pybedtools
import os

# import glob
# import shutil
original_species = "Mm"
wdir = os.getcwd() + "/"
protein = "Cenpa"


# Hard code "reference_contigs_dict" for human and mouse, potentially other main species
# Cose you need the user to have it with the packcage without digging to instal folder
# Done

# Provide the user with a list of possible genes per assembly (to avoid wrong naming)
# Check protein list for valid gene names
# Write number of transcript into the cds header
# Write a test function that would ensure that all first exons start with ATG and all the last exons end with TGA, TAG and TAA
# Look at Cenpa, Cenpb, Cenpc, Cenpe, Cenpf in caroli and pahari and check with ensembl cds for identity

# CENP-A has UTR_3 spanning last two exons -> does not work with the current extract_exons function!

rules = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
         "Y": "R", "R": "Y", "W": "W", "S": "S", "K": "M", "M": "K",
         "D": "H", "H": "D", "V": "B", "B": "V", "X": "X"}

def generate_reference_genome_object(wdir, original_species):
    """Generates a reference Genome object using pyensembl as a wrapper for ensembl database"""
    
    # THE ORIGINAL SPECIES NAME SHOULD BE STREAMLINED TO "Hs" etc...
    
    # default is mouse
    mouse_names = {"Mm", "mouse", "Mus musculus", "mus musculus"}
    human_names = {"Hs", "human", "Homo sapiens", "homo sapiens"}
    if original_species in mouse_names:
        #header = "musculus"
        species = "mus musculus"
        release = 100
        reference_genome_contigs_dict = reference_genome_dict_generator.get_reference_genome_contigs_dict(original_species)
        reference_genome_name = "MUSCULUS_genome"
    
    elif original_species in human_names:
        #header = "sapiens"
        species = "homo sapiens"
        release = 100
        reference_genome_contigs_dict = reference_genome_dict_generator.get_reference_genome_contigs_dict(original_species)
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
    
    # match contigs named between pyensembl Genome object and GeneBank assembly
    #reference_contigs_dict = {}
    #with open(reference_genome_contigs_file, "r") as f:
    #    file = f.readlines()
    #    for line in file:
    #        l = line.split("\t")
    #        reference_contigs_dict[l[0]] = l[4]
    
    return reference_genome_present, ensembl, original_species, reference_genomes_path, \
            reference_genome_name, reference_genome_contigs_dict, biotype


def extract_input(wdir, original_species, reference_genome_name, reference_genomes_path, 
                  reference_genome_contigs_dict, ensembl, biotype, protein):
    """Extracts all input sequences required by FREEDA from the indicated reference genome (original_species)"""
    
    
    blast_input_path = wdir + "Blast_input/"
    coding_sequence_input_path = wdir + "Coding_sequences/"
    gene_input_path = wdir + "Genes/"
    exons_input_path = wdir + "Exons/"

    # find coding sequence
    transcript, longest_transcript_id, gene_id, contig, strand, UTR_5, UTR_3, cds_sequence = extract_cds(ensembl, 
                original_species, coding_sequence_input_path, protein, biotype)
    # find protein sequence
    extract_protein(original_species, blast_input_path, protein, strand, transcript)
    # find gene sequence
    extract_gene(original_species, gene_input_path, ensembl, contig, strand, gene_id, 
                reference_genomes_path, reference_genome_name, reference_genome_contigs_dict, protein)
    # find exons sequence
    extract_exons(original_species, protein, exons_input_path, reference_genomes_path, reference_genome_name, 
                  contig, strand, transcript, reference_genome_contigs_dict, UTR_5, UTR_3, cds_sequence)
    
    
def extract_protein(original_species, blast_input_path, protein, strand, transcript):
    """Extracts protein sequence based on Transcript object"""
    
    # find proteins sequence
    protein_sequence = transcript.protein_sequence
    # get and save the sequence
    parse_sequence(original_species, blast_input_path, protein_sequence, 
                   protein, strand = None, sequence_type = "protein")
    
    
def extract_exons(original_species, protein, exons_input_path, reference_genomes_path, reference_genome_name, 
                  contig, strand, transcript, reference_genome_contigs_dict, UTR_5, UTR_3, cds_sequence):
    """Extracts exoms sequence based on Transcript object"""
    
    # THIS FUNCTION WAS NOT TESTED FOR UTR_3 SPANNING MORE THAN ONE EXON!
        
    # get start and stop for each exon
    exon_intervals = transcript.exon_intervals
   
    # get sequence of each exon
    exons_filename = protein + "_" + original_species + "_exons.fasta"
    exons_file = open(exons_input_path + exons_filename, "w")
    
    number = 0
    last_exon_number = 1000
    single_coding_exon = False
    coding_exons = {}

    exon_intervals_copy = exon_intervals[::]
    
    # Get first coding exon
    for exon in exon_intervals_copy:
        exon_fasta_sequence = get_single_exon(original_species, protein, reference_genome_contigs_dict, 
            reference_genomes_path, reference_genome_name, rules, contig, strand, number, exon, last_exon_number)
        
        if len(UTR_5) > 0:
            
            # This is a non-coding exon
            if len(UTR_5) >= len(exon_fasta_sequence):
                # trim the UTR_5 and dont write this exon into exons file
                UTR_5 = UTR_5[len(exon_fasta_sequence):]
                # print("Trimmed UTR_5: %s" % UTR_5)
                # reset exon count
                number == 0
                # remove that exon from the intervals list
                exon_intervals = exon_intervals[1:]
                continue
            
            # This is the first coding exon
            else:
                number += 1
                # trim and write exon -> already in correct orientation
                exon_fasta_sequence_trimmed = exon_fasta_sequence.replace(UTR_5, "")
                exon_fasta_sequence = exon_fasta_sequence_trimmed
                UTR_5 = 0
                # remove that exon from the intervals list if there are more than 1 exons listed
                if len(exon_intervals) != 1:
                    # write header into coding exons dict (number should = 1)
                    coding_exons[number] = exon_fasta_sequence
                    exon_intervals = exon_intervals[1:]
                else:
                    # This is NOT triggred if UTR_3 is divided into more than one exon!
                    single_coding_exon = True
                
                break
    
    # Get last coding exon
    for exon in reversed(exon_intervals_copy):
        
        # If there is more than one coding exon
        if single_coding_exon == False:
            last_exon_number -= 1
            exon_fasta_sequence = get_single_exon(original_species, protein, reference_genome_contigs_dict, 
                reference_genomes_path, reference_genome_name, rules, contig, strand, number, exon, last_exon_number)
        
            if len(UTR_3) > 0:
            
                # This is a non-coding exon
                if len(UTR_3) > len(exon_fasta_sequence):
                    # trim the UTR_3 and dont write this exon into exons file
                    UTR_3 = UTR_3[:len(exon_fasta_sequence)]
                    # remove that exon from the intervals list
                    exon_intervals = exon_intervals[:-1]
                    continue
            
                # This is the last coding exon
                else:
                    # trim and write exon -> already in correct orientation
                    exon_fasta_sequence_trimmed = exon_fasta_sequence.replace(UTR_3, "")
                    exon_fasta_sequence = exon_fasta_sequence_trimmed
                    UTR_3 = 0
                    # write header into coding exons dict -> pick a huge arbitrary placeholder number
                    coding_exons[float("inf")] = exon_fasta_sequence
                    # remove that exon from the intervals list
                    exon_intervals = exon_intervals[:-1]
                    break
        
        # If there is only one coding exon
        else:
            number == 1
            # trim and write exon sequence coming from UTR_5 trimming -> already in correct orientation
            exon_fasta_sequence_trimmed = exon_fasta_sequence.replace(UTR_3, "")
            exon_fasta_sequence = exon_fasta_sequence_trimmed
            UTR_3 = 0
            # write header into coding exons dict (number should = 1)
            coding_exons[number] = exon_fasta_sequence
            break
            
               
    # Get the rest of the coding exons
    if single_coding_exon == False:
        for exon in exon_intervals:
            exon_fasta_sequence = get_single_exon(original_species, protein, reference_genome_contigs_dict, 
                reference_genomes_path, reference_genome_name, rules, contig, strand, number, exon, last_exon_number)
        
            # write header into coding exons dict
            number += 1
            coding_exons[number] = exon_fasta_sequence
        
    # Reset the number of the last exon
    final_order = sorted([number for number, exon in coding_exons.items()])
    final_coding_exons = {number : coding_exons[number] for number in final_order}
    
    #print(coding_exons)
    #print(final_order)
    #print(final_coding_exons)
    
    for nr, exon_sequence in final_coding_exons.items():
        if nr != float("inf"):
            header = ">" + protein + "_exon_" + str(nr)
        else:
            header = ">" + protein + "_exon_" + str(len(final_coding_exons))
        exons_file.write(header + "\n")
        exons_file.write(exon_sequence + "\n")
    
    exons_file.close()
            

def get_single_exon(original_species, protein, reference_genome_contigs_dict, 
        reference_genomes_path, reference_genome_name, rules, contig, strand, number, exon, last_exon_number):
    """Returns a header and a string representation of a current exon."""
    
    if last_exon_number < 1000:
        number = last_exon_number
    else:
        number += 1
        
    # need to provide an off by 1 offset for start on both strands
    start = exon[0] - 1
    stop = exon[1]
    header = protein + "_" + original_species + "_exon_" + str(number)
    exon_bed_filename = header + ".bed"
    exon_fasta_filename = header + ".fasta"
    with open(exon_bed_filename, "w") as b:
        b.write(reference_genome_contigs_dict[contig] + "\t")
        b.write(str(start) + "\t")
        b.write(str(stop))
        
    bed_object = pybedtools.BedTool(exon_bed_filename)
    bed_object = bed_object.sequence(fi = reference_genomes_path + \
                        reference_genome_name + ".fasta")
    # makes a file and writes a fasta seq into it
    exon_fasta_sequence = bed_object.save_seqs(exon_fasta_filename)
    
    with open(exon_fasta_filename, "r") as f:
        file = f.readlines()
        exon_fasta_sequence = file[1].rstrip("\n").upper()
    
    if strand == "-":
        complement = ""
        for base in exon_fasta_sequence:
            complement = rules[base] + complement
        exon_fasta_sequence = complement
    
    os.remove(exon_bed_filename)
    os.remove(exon_fasta_filename)
    
    return exon_fasta_sequence

def extract_gene(original_species, gene_input_path, ensembl, contig, strand, gene_id, 
        reference_genomes_path, reference_genome_name, reference_genome_contigs_dict, protein):
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
        
    matched_contig = reference_genome_contigs_dict[str(contig)]
    
    # make a bedtool file
    with open(gene_bed_filename, "w") as b:
        b.write(str(matched_contig) + "\t" + str(start) + "\t" + str(end))
    
    # parse that bedtool file using BedTools object
    bed_object = pybedtools.BedTool(gene_bed_filename)
    gene_fasta_sequence = bed_object.sequence(fi = reference_genomes_path + \
                            reference_genome_name + ".fasta")
        
    os.remove(gene_bed_filename)
        
    parse_sequence(original_species, gene_input_path, gene_fasta_sequence, 
                   protein, strand, sequence_type="gene")


def parse_sequence(original_species, output_path, fasta_sequence, protein, strand, sequence_type):
    """Parses either a pybedtools.bedtool.BedTool object or string into a file"""
    
    # check if its a bedtool object
    if isinstance(fasta_sequence, pybedtools.bedtool.BedTool):
        file = open(fasta_sequence.seqfn)
        content = file.readlines()
        header = ">" + protein + "_" + original_species + "_" + sequence_type
        sequence = content[1].upper().rstrip("\n")
        file.close()
    
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
    parse_sequence(original_species, coding_sequence_input_path, cds_sequence, 
                   protein, strand, sequence_type="cds")
    
    if start_codon_present and stop_codon_present:
        return transcript, longest_transcript_id, gene_id, contig, strand, UTR_5, UTR_3, cds_sequence
    else:
        print("Chosen transcript for protein: %s does not have START and STOP annotated" % protein)
        

    