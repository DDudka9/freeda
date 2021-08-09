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
from bioservices import UniProt
from Bio import SeqIO
import pyensembl
import pybedtools
import os
import glob
import subprocess
import shutil
import requests
import tarfile

# import glob
# import shutil
original_species = "Mm"
wdir = os.getcwd() + "/"
protein = "Haus1"
#reference_genome_name = "MUSCULUS_genome"


# Hard code "reference_contigs_dict" for human and mouse, potentially other main species
# Cose you need the user to have it with the packcage without digging to instal folder
# Done

# Check protein list for valid gene names
# DONE

# Write number of transcript into the cds header
# DONE

# Write a test function that would ensure that all first exons start with ATG and all the last exons end with TGA, TAG and TAA
# DONE

# Look at Cenpa, Cenpb, Cenpc, Cenpe, Cenpf in caroli and pahari and check with ensembl cds for identity

# Make a function that aligns protein sequence from pyensembl and from AlphaFold

# Make a function deleting microexon from exons
# DONE -> remove microexons check from further functions? -> exon finder doesnt know one exon was skipped!!!
# Still there is one additional exon found.... cose exons are first aligned using CDS which still contains microexon
# FURTHER CHECKPOINT MIGHT BE TRIGGERED IF FINAL CDS IS OUT OF FRAME (when translated???)

# Think how to deal with gene_and_cds_reader module -> also detects microexons etc -> leave it in case of manual input
# but no microexons would be detected and no warnings issued (which is a problem)

# Also stucture model should be input together with the rest of the input

rules = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
         "Y": "R", "R": "Y", "W": "W", "S": "S", "K": "M", "M": "K",
         "D": "H", "H": "D", "V": "B", "B": "V", "X": "X"}


def get_uniprot_id(wdir, original_species, protein):
    """Retrieves all possible uniprot ids to be matched against structure prediction from AlphaFold"""
    
    if original_species in {"Mm", "Mouse", "mouse", "Mus musculus", "mus musculus"}:
        original_species_number = "10090"
    
    if original_species in {"Hs", "Human", "human", "Homo sapiens", "homo sapiens"}:
        original_species_number = "9606"
    
    possible_uniprot_ids = set()
    
    # generate a search for given protein
    u = UniProt(verbose=False)
    print("\nExtracting valid uniprot ids for protein: %s ...\n" % protein)
    data = u.search(protein + "+and+taxonomy:" + original_species_number, frmt="tab", limit=5,
             columns="genes,length,id")
    data_list = data.split("\n")
    for line in data_list:
        elements = line.split("\t")
        if protein.upper() in elements[0].upper():
            possible_uniprot_ids.add(elements[2])
    
    return possible_uniprot_ids


def fetch_structure_prediction(wdir, original_species, protein, possible_uniprot_ids):
    """Finds and extracts structure prediction model from AlphaFold database for the species"""
    
    if original_species in {"Mm", "Mouse", "mouse", "Mus musculus", "mus musculus"}:
        handle = "MOUSE"
    if original_species in {"Hs", "Human", "human", "Homo sapiens", "homo sapiens"}:
        handle = "HUMAN"
        
    # make folder to host the structure prediction
    structure_path = wdir + "Structures/" + protein + "_" + original_species
    if not os.path.isdir(structure_path):
        os.makedirs(structure_path)
    
    # find prediction model based on possible uniprot ids for the protein
    path_to_original_species_structures = [path for path in glob.glob(wdir + "*") 
                                                   if path.endswith(handle + ".tar")][0]
    tar = tarfile.open(path_to_original_species_structures)
    print("Looking for structure in AlphaFold database for: %s ...\n" % protein)
    all_structures = [name for name in tar.getnames() if name.endswith(".pdb.gz")]
    # this assumes that there is only one predicion model for a given protein in AlpghaFold
    all_structures_compressed = [s for s in all_structures if s.split("-")[1] in possible_uniprot_ids]

    
    # try finding the structure
    try:
        structure_filename_compressed = all_structures_compressed[0]
        uniprot_id = structure_filename_compressed.split("-")[1]
        
        # extract a compressed model
        path_to_tarfile = wdir + structure_filename_compressed
        tar.extract(structure_filename_compressed)
        
        # need to use subprocess to decompress (tarfile or gunzip modules dont work directly)
        cmd = ["gunzip", path_to_tarfile] # deletes zipped file automatically
        subprocess.call(cmd)
        
        # define path to decompressed file and move to Structures folder
        structure_filename_decompressed = structure_filename_compressed.replace(".gz", "")
        path_to_tarfile_decompressed = wdir + structure_filename_decompressed
        
        # silence warnings from Biopython about missing header in model (has to follow import)
        import warnings
        from Bio import BiopythonWarning
        warnings.simplefilter('ignore', BiopythonWarning)
        
        # get model sequence and length
        with open(path_to_tarfile_decompressed, 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                model_seq = record.seq
        
        print("Based on uniprot id: %s structure prediction for protein: %s has been found (%s aa)" 
                                                  % (uniprot_id, protein, len(model_seq)))
        shutil.move(path_to_tarfile_decompressed, structure_path + "/" + structure_filename_decompressed)
        
        return model_seq
    
    # didnt find structure
    except IndexError:
        print("...WARNING...: Structure prediction for protein: %s HAS NOT BEEN FOUND -> cannot run PyMOL\n" % protein)
        print("...SUGGESTION...: You can use your own model (ex. from PDB; fragments are ok but no sequence mismatches!)\n")
        
        return False

# THIS IS NOT NEEDED ANYMORE
def get_prediction(wdir, original_species, protein):
    """Generates an Alpha Fold url for a given protein based on its UniProt id"""
    
    # check if the protein structure folder exists
    structure_dir = wdir + "Structures/" + protein + "_" + original_species
    if os.path.isdir(structure_dir):
        # if structure prediction already exists - skip it
        if len(os.listdir(structure_dir)) != 0:
            return True
        else:
            pass
    else:  
        os.makedirs(wdir + "Structures/" + protein + "_" + original_species)
    
    supported_species_names = {"Mm" : "mouse", "Hs" : "human"}
    organism = supported_species_names.get(original_species) 
    
    url = "https://www.uniprot.org/uniprot/?query=" + protein + "+organism:" + organism + "&sort=score&columns=id,reviewed,genes,organism&format=tab" 
    downloaded_obj = requests.get(url)
    filename = protein + ".txt"

    with open(filename, "wb") as file:
        file.write(downloaded_obj.content)
    
    with open(filename, "r") as file:
        f = file.readlines()
        for line in f:
            columns = line.split("\t")
            list_of_names = columns[2].split(" ")
            # in case there are multiple gene names
            if len(list_of_names) > 1 and protein in list_of_names:
                uniprot_id = line.split("\t")[0]
                break
            else:
                # take the most likely id (reviewed are favored)
                if protein in list_of_names:
                    uniprot_id = columns[0]
                    break

    prediction_url = "https://www.alphafold.ebi.ac.uk/entry/" + str(uniprot_id)
        
    os.remove(filename)
    
    return prediction_url
        

def generate_reference_genome_object(wdir, original_species, reference_genome_name):
    """Generates a reference Genome object using pyensembl as a wrapper for ensembl database"""

    mouse_names = {"Mm", "Mouse", "mouse", "Mus musculus", "mus musculus"}
    human_names = {"Hs", "Human", "human", "Homo sapiens", "homo sapiens"}
    
    if original_species in human_names:
        species = "homo sapiens"
        release = 100
        reference_genome_contigs_dict = reference_genome_dict_generator.get_reference_genome_contigs_dict(original_species)
    
    # default is mouse for now
    else:
        species = "mus musculus"
        release = 100
        reference_genome_contigs_dict = reference_genome_dict_generator.get_reference_genome_contigs_dict(original_species)
    
    
    # make sure original species genome (reference genome) is present
    reference_genomes_path = wdir + "Reference_genomes/"
    # check if reference genome is present -> exit 
    reference_genome_present = check_genome_present(reference_genomes_path, reference_genome_name, reference_genome=True)

    # get assembly database
    ensembl = pyensembl.EnsemblRelease(release, species)
    # define biotype (FREEDA deals only with protein coding sequences)
    biotype = "Protein coding"
    
    return reference_genome_present, ensembl, original_species, reference_genomes_path, \
            reference_genome_contigs_dict, biotype


def extract_input(wdir, original_species, reference_genome_name, reference_genomes_path, 
                  reference_genome_contigs_dict, ensembl, biotype, protein, model_seq):
    """Extracts all input sequences required by FREEDA from the indicated reference genome (original_species)"""
    
    model_matches_input = False
    microexon_present = False
    
    blast_input_path = wdir + "Blast_input/"
    coding_sequence_input_path = wdir + "Coding_sequences/"
    gene_input_path = wdir + "Genes/"
    exons_input_path = wdir + "Exons/"

    # check if provided protein.txt list has valid gene names
    all_genes = {g for g in ensembl.gene_names()}
    if protein not in all_genes:
        input_correct = False
        print("\n...WARNING...: Reference genome lacks gene name: "'%s'"\n" % protein)
        return input_correct, model_matches_input, microexon_present
    
    # find coding sequence
    transcript, selected_transcript_id, gene_id, contig, strand, UTR_5, UTR_3, cds_sequence_expected, matching_length = extract_cds(ensembl, 
                original_species, coding_sequence_input_path, protein, biotype, model_seq)
    # find protein sequence
    model_matches_input = extract_protein(original_species, blast_input_path, protein, strand, transcript, model_seq, matching_length)
    # find gene sequence
    extract_gene(original_species, gene_input_path, ensembl, contig, strand, gene_id, 
                reference_genomes_path, reference_genome_name, reference_genome_contigs_dict, protein, transcript)
    # find exons sequence
    input_correct, microexon_present = extract_exons(wdir, original_species, protein, exons_input_path, reference_genomes_path, reference_genome_name, 
                  contig, strand, transcript, reference_genome_contigs_dict, UTR_5, UTR_3, cds_sequence_expected)
    
    return input_correct, model_matches_input, microexon_present
    

def extract_protein(original_species, blast_input_path, protein, strand, transcript, model_seq, matching_length):
    """Extracts protein sequence based on Transcript object"""
    
    model_matches_input = False
    
    # find proteins sequence
    protein_sequence = transcript.protein_sequence
    # get and save the sequence
    parse_sequence(original_species, blast_input_path, protein_sequence, 
                   protein, transcript, strand = None, sequence_type = "protein")
    # check if protein length and sequence match that of the model
    if matching_length == True:
        if model_seq == protein_sequence:
            model_matches_input = True
            
    return model_matches_input
    
    
def extract_exons(wdir, original_species, protein, exons_input_path, reference_genomes_path, reference_genome_name, 
                  contig, strand, transcript, reference_genome_contigs_dict, UTR_5, UTR_3, cds_sequence_expected):
    """Extracts exoms sequence based on Transcript object"""
    
    # get start and stop for each exon
    exon_intervals = transcript.exon_intervals

    start_codon_offset = transcript.first_start_codon_spliced_offset
   
    # get sequence of each exon
    exons_filename = protein + "_" + original_species + "_exons.fasta"
    exons_file = open(exons_input_path + exons_filename, "w")
    
    number = 0
    coding_exons = {}
    UTR_5_length = len(UTR_5)
    UTR_3_length = len(UTR_3)
    removed_on_coding_from_start = False
    removed_on_coding_from_end = False
    start_codon_present = False
    stop_codon_present = False
    nr_of_removed_non_coding_exons_start = 0
    nr_of_removed_non_coding_exons_end = 0
    cds_from_exons = ""
    microexon_present = False
    input_correct = True
    
    # make a true copy of exon coordinates
    exon_intervals_copy = exon_intervals[::]
    
    # Get first coding exon
    for exon in exon_intervals_copy:
        
        if exon[1]-exon[0] > UTR_5_length:
            #print("start_offset_final: %s\n" % start_codon_offset)
            break
        
        # trim off the non-coding exons containing UTR_5
        if exon[1]-exon[0] <= UTR_5_length:
            removed_on_coding_from_start = True
            nr_of_removed_non_coding_exons_start += 1
            exon_intervals = exon_intervals[1:]
            trimmed_UTR_5_length = UTR_5_length - (exon[1]-exon[0])
            start_codon_offset = trimmed_UTR_5_length
            UTR_5_length = trimmed_UTR_5_length
    
    # Get last coding exon
    for exon in reversed(exon_intervals_copy):
        
        if exon[1]-exon[0] > UTR_3_length:
            break
        
        # trim off the non-coding exons containing UTR_3
        if exon[1]-exon[0] <= UTR_3_length:
            removed_on_coding_from_end = True
            nr_of_removed_non_coding_exons_end += 1
            exon_intervals = exon_intervals[:-1]
            trimmed_UTR_3_length = UTR_3_length - (exon[1]-exon[0])
            UTR_3_length = trimmed_UTR_3_length
    
    # Get all coding exons
    for exon in exon_intervals:
        number += 1
        start = exon[0] - 1
        stop = exon[1]
        
        if number == 1 and strand == "+":
            # WORKS OK? (tested on Cenpa)
            if removed_on_coding_from_start == False:
                start = start + start_codon_offset
            # WORKS OK? (tested on Cenpi)
            else:
                start = start + start_codon_offset - nr_of_removed_non_coding_exons_start
        
        if number == 1 and strand == "-":
            # WORKS OK? (tested on Cenpr -> Itgb3bp)
            if removed_on_coding_from_start == False:
                stop = stop - start_codon_offset
            # WORKS OK? (tested on Cenpo)
            else:
                stop = stop - start_codon_offset + nr_of_removed_non_coding_exons_start
        
        if number == len(exon_intervals) and strand == "+":
            # WORKS OK? (tested on Cenpi)
            if removed_on_coding_from_end == False:
                stop = stop - UTR_3_length
            # WORKS OK? (tested on Cenpa)
            else:
                stop = stop - UTR_3_length + nr_of_removed_non_coding_exons_end
        
        if number == len(exon_intervals) and strand == "-":
            # WORKS OK? (tested on Cenpo)
            if removed_on_coding_from_end == False:
                start = start + UTR_3_length
            # WORKS OK? (tested on Cenpr -> Itgb3bp)
            else:
                start = start + UTR_3_length - 1

        exon_fasta_sequence = get_single_exon(original_species, protein, reference_genome_contigs_dict, 
                reference_genomes_path, reference_genome_name, rules, contig, strand, number, exon, start, stop)
        coding_exons[number] = exon_fasta_sequence

    for nr, exon_sequence in coding_exons.items():
        
        if len(exon_sequence) < 20:
            microexon_present = True
            print("\n...WARNING...:Exon %s in: %s is a microexon (%sbp; hard to align) " \
                    "-> eliminated from exons and cds\n" % (str(nr), protein, len(exon_sequence)))
            # check if the last exon (microexon) had a stop codon before its trimmed (ex. Pot1a)
            if nr == 1 and exon_sequence.startswith("ATG"):
                start_codon_present = True
            # check if the last exon (microexon) had a stop codon before its trimmed (ex. Numa1)
            if nr == len(coding_exons) and exon_sequence[-3:] in ["TGA", "TAG", "TAA"]:
                stop_codon_present = True
            continue
            
        else:
            cds_from_exons += exon_sequence
        
            if nr == 1 and exon_sequence.startswith("ATG"):
                start_codon_present = True
            if nr == len(coding_exons) and exon_sequence[-3:] in ["TGA", "TAG", "TAA"]:
                stop_codon_present = True
        
        header = ">" + transcript.name + "_" + original_species + "_exon_" + str(nr)
        exons_file.write(header + "\n")
        exons_file.write(exon_sequence + "\n")
        
        
    if start_codon_present == False:
        print("\nFirst exon in: %s in missing a START codon!!!\n" % header.split("_")[0].replace(">",""))
        input_correct = False
        
    if stop_codon_present == False:
        print("\nLast exon in: %s in missing a STOP codon!!!\n" % header.split("_")[0].replace(">",""))
        input_correct = False
    
    # use the coding sequence assembled from exons to overwrite the cds_expected_sequence
    if microexon_present == True:
        output_path = wdir + "Coding_sequences/"
        parse_sequence(original_species, output_path, cds_from_exons, protein, transcript, strand, sequence_type="cds")

    if microexon_present == False and cds_from_exons != cds_sequence_expected:
        print("\nExons FAILED to assemble expected CDS for: %s\n" % header.split("_")[0].replace(">",""))
        input_correct = False
        
    exons_file.close()

    return input_correct, microexon_present

def get_single_exon(original_species, protein, reference_genome_contigs_dict, reference_genomes_path, 
                    reference_genome_name, rules, contig, strand, number, exon, start, stop):
    """Returns a header and a string representation of a current exon."""

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
        reference_genomes_path, reference_genome_name, reference_genome_contigs_dict, protein, transcript):
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
                   protein, transcript, strand, sequence_type="gene")


def parse_sequence(original_species, output_path, fasta_sequence, protein, transcript, strand, sequence_type):
    """Parses either a pybedtools.bedtool.BedTool object or string into a file"""
    
    # check if its a bedtool object
    if isinstance(fasta_sequence, pybedtools.bedtool.BedTool):
        file = open(fasta_sequence.seqfn)
        content = file.readlines()
        header = ">" + transcript.name + "_" + original_species + "_" + sequence_type
        sequence = content[1].upper().rstrip("\n")
        file.close()
    
    # or its a plain string
    if not isinstance(fasta_sequence, pybedtools.bedtool.BedTool):
        header = ">" + transcript.name + "_" + original_species + "_" + sequence_type
        sequence = fasta_sequence
    
    # or plain sting on reverse string (gene only)
    if sequence_type == "gene" and strand == "-":
        complement = ""
        for base in sequence:
            complement = rules[base] + complement
        sequence = complement
    
    # write it into a file
    filename = protein + "_" + original_species + "_" + sequence_type
    with open(output_path + filename + ".fasta", "w") as f:
        f.write(header + "\n")
        f.write(sequence)


def extract_cds(ensembl, original_species, coding_sequence_input_path, protein, biotype, model_seq):
    """Extracts coding sequence by creating Transcript object"""
    
    all_transcripts_ids = ensembl.transcript_ids_of_gene_name(protein)
    all_transcripts_dict = {}

    # find all possible transcripts
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
            cds_sequence_expected = transcript.coding_sequence
        except ValueError:
            cds_sequence_expected = ""
            all_transcripts_dict[t].append(cds_sequence_expected)
        finally:
            all_transcripts_dict[t].append(cds_sequence_expected)
        
        # get length of cds
        length = len(cds_sequence_expected)
        all_transcripts_dict[t].append(length)
    
    # ask user if they have any preference
    #user_dict = {}
    #for t, features in all_transcripts_dict.items():
    #    length = (features[-1]-3) * 3
    #    if length > 0:
    #        user_dict[features[3]] = str(int((features[-1]-3) / 3)) + " aa"
    
    #tries = 3
    #preference = ""
    #while tries > 0:
    #    tries -= 1
    #    preference = str(input("Choose best transcript for your analysis (strongly recommended) or press ENTER (longest cds):\n %s\n" 
    #                       % user_dict))
    #    if preference in user_dict or preference == "":
    #        tries = 0
    
    # pick the trancript of preference (recommended)
    #if preference in user_dict:
    #    for t, features in user_dict.items():
    #        if t == preference:
    #            for t, features in all_transcripts_dict.items():
    #                if features[3] == preference:
    #                    longest_transcript_id = t
    #                    length = features[-1]
    
    # compare length of translated coding sequence with model amino acid seuqence length
    matching_length = False
    for t, features in all_transcripts_dict.items():
        cds_length_aa = (features[-1]-3)/3
        if len(model_seq) == cds_length_aa:
            selected_transcript_id = t
            matching_length = True
    
    # if there is no transcript of preference, pick the one with longest cds (not recommended)
    if matching_length == False:
        selected_transcript_id = None
        length = 0
        for t, features in all_transcripts_dict.items():
            if features[-1] > length:
                selected_transcript_id = t
                length = features[-1]
    
    # unpack features of the longest transcript
    gene_id, contig, strand, transcript_name, start, end, cds_sequence_expected, length = all_transcripts_dict[selected_transcript_id]
    # rebuild Transcript object based on the longest transcript
    transcript = pyensembl.Transcript(selected_transcript_id, transcript_name, 
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
    parse_sequence(original_species, coding_sequence_input_path, cds_sequence_expected, 
                   protein, transcript, strand, sequence_type="cds")
    
    if start_codon_present and stop_codon_present:
        return transcript, selected_transcript_id, gene_id, contig, strand, UTR_5, UTR_3, cds_sequence_expected, matching_length
    else:
        print("Chosen transcript for protein: %s does not have START and STOP annotated" % protein)
        
