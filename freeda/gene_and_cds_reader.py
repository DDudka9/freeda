#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:28:24 2021

@author: damian

Finds and parses a cds and gene sequences for a given protein from the original species

"""

import re
import logging

def find_gene_and_cds(wdir, protein_name, original_species):
    
    Mm_exons, expected_exons, microexons, microexons_seqs = get_Mm_exons(wdir, protein_name, original_species)

    # open according cds fasta file
    with open(wdir + "Coding_sequences/" + protein_name + "_" + original_species + "_cds.fasta", "r") as f:
        make_linear = [line.rstrip("\n") for line in f.readlines()]
        for line in make_linear:
            if line.startswith(">"):
                cds_linear = line + "\n"
            else:
                cds_linear += line
    
    cds = cds_linear
                
    if microexons != []:
        # grab sequence of each microexom
        to_delete = [seq for seq in microexons_seqs]
            
        # remove these sequenses from the coding sequence
        for seq in to_delete:
            cds_edited = cds_linear.replace(seq, "")
            cds = cds_edited

    # open according gene fasta file
    with open(wdir + "Genes/" + protein_name + "_" + original_species + "_gene.fasta", "r") as f:
        gene = f.read()
        
    return cds, gene, Mm_exons, microexons, expected_exons


def get_Mm_exons(wdir, protein_name, original_species): # works well -> use for cloning after synteny check
    # get path to the exons for given protein
    
    Mm_exons = {}
    count = 0
    seq_recorded = False
    header = ""
    seq = ""
    microexons = []
    
    with open(wdir + "Exons/" + protein_name + "_" + original_species + "_exons.fasta", "r") as f:
        file = f.read()

        for line in re.split("\n", file):
            
            # this statement executes last
            if line.startswith(">") and seq_recorded == True:
                # record header, sequence and length of the exon
                Mm_exons[count] = (header, seq, len(seq))
                seq_recorded = False
                seq = ""
            
            # this statement executes first
            if line.startswith(">") and seq_recorded == False:
                count += 1
                head = line.lstrip(">").rstrip("\n")
                header = ">_" + "exon_" + str(count) + "_" + head
                seq_recorded = True
            
            # this statement executes next
            if not line.startswith(">"):
                seq = seq + line.rstrip("\n").upper()
            
        # record the last exon
        Mm_exons[count] = (header, seq, len(seq))
        
        # double check if all exons together are in frame (they should be)
        exon_total_length = 0
        for number, v  in Mm_exons.items():
            exon_total_length += v[2]
            
        if exon_total_length % 3 != 0:
            message = "\nCDS of %s in original species is NOT in frame."\
                % protein_name
            print(message)
            logging.info(message)
         
        # detect possible microexons (hardcoded 20bp limit) -> too difficult to align
        
        # CHANGED LIMIT TO 18bp 06/05/2021 to run CENPC exon 1
        microexons = []
        microexons_seqs = []
        minimum_exon_length = 18
        editing_threshold = 9
        for exon in Mm_exons:
            length = len(Mm_exons[exon][1])
            if length < minimum_exon_length:
                microexons_seqs.append(Mm_exons[exon][1])
                if exon != list(Mm_exons)[-1]:
                    microexons.append(exon)
                    message = "**** Exon nr " + str(exon) + " is a microexon: " + str(length) + \
                        " bp (min. length = %sbp)" % minimum_exon_length
                    print(message)
                    logging.info(message)
                
                else:
                    message = "**** Last exon nr " + str(exon) + " is a microexon: " + str(length) + \
                        " bp -> allowed in last exons (min. length for other exons = %sbp)" % minimum_exon_length
                    print(message)
                    logging.info(message)
                    
            if length < editing_threshold:
                message = "**** CAUTION - Exon nr " + str(exon) + " is a VERY SHORT microexon: " + str(length) + \
                    " bp (please manually remove it from the original species CDS; automatic editing advised >= %sbp)" % editing_threshold
                print(message)
                logging.info(message)
                
        # eliminate microexons
        if microexons != []:
            for microexon in microexons:
                del Mm_exons[microexon]
                message = "**** Exon nr " + str(microexon) + " was skipped to ease alignment"
                print(message)
                logging.info(message)
            # re-number exons (enumerate method starting from 1)
            #M = {}
            #for number, features in enumerate(Mm_exons, 1):
            #    M[number] = Mm_exons[features]
            #Mm_exons = M
        
        expected_exons = tuple(e for e, features in Mm_exons.items())
        message = "........................................................................\n\n" \
                "ANALYZING PROTEIN: %s \n\n" \
                "........................................................................\n\n" \
                "Expected exons : %s" % (protein_name, str(expected_exons))
        print(message)
        logging.info(message)
        
    return Mm_exons, expected_exons, microexons, microexons_seqs



"""

def get_Mm_exons(wdir, protein_name): # works well -> use for cloning after synteny check
    # get path to the exons for given protein
    files = []
    files.append(glob.glob(wdir + "Exons/" + "*" + protein_name + "_*.fasta")[0])

    # sort exon names and sequences into a dictionary
    for f in files:
        file = open(f, "r").read()
        Mm_exons = {}
        count = 0
        seq_recorded = False
        header = ""
        seq = ""
        microexons = []

        for line in re.split("\n", file):
            
            # this statement executes last
            if line.startswith(">") and seq_recorded == True:
                # record header, sequence and length of the exon
                Mm_exons[count] = (header, seq, len(seq))
                seq_recorded = False
                seq = ""
            
            # this statement executes first
            if line.startswith(">") and seq_recorded == False:
                count += 1
                head = line.lstrip(">").rstrip("\n")
                header = ">_" + "exon_" + str(count) + "_" + head
                seq_recorded = True
            
            # this statement executes next
            if not line.startswith(">"):
                seq = seq + line.rstrip("\n").upper()
            
        # record the last exon
        Mm_exons[count] = (header, seq, len(seq))
        
        # double check if all exons together are in frame (they should be)
        exon_total_length = 0
        for number, v  in Mm_exons.items():
            exon_total_length += v[2]
            
        if exon_total_length % 3 != 0:
            message = "\nCDS of %s in original species is NOT in frame."\
                % protein_name
            print(message)
            logging.info(message)
         
        # detect possible microexons (hardcoded 20bp limit) -> too difficult to align
        
        # CHANGED LIMIT TO 18bp 06/05/2021 to run CENPC exon 1
        microexons = []
        microexons_seqs = []
        minimum_exon_length = 18
        editing_threshold = 9
        for exon in Mm_exons:
            length = len(Mm_exons[exon][1])
            if length < minimum_exon_length:
                microexons_seqs.append(Mm_exons[exon][1])
                if exon != list(Mm_exons)[-1]:
                    microexons.append(exon)
                    message = "**** Exon nr " + str(exon) + " is a microexon: " + str(length) + \
                        " bp (min. length = %sbp)" % minimum_exon_length
                    print(message)
                    logging.info(message)
                
                else:
                    message = "**** Last exon nr " + str(exon) + " is a microexon: " + str(length) + \
                        " bp -> allowed in last exons (min. length for other exons = %sbp)" % minimum_exon_length
                    print(message)
                    logging.info(message)
                    
            if length < editing_threshold:
                message = "**** CAUTION - Exon nr " + str(exon) + " is a VERY SHORT microexon: " + str(length) + \
                    " bp (please manually remove it from the original species CDS; automatic editing advised >= %sbp)" % editing_threshold
                print(message)
                logging.info(message)
                
        # eliminate microexons
        if microexons != []:
            for microexon in microexons:
                del Mm_exons[microexon]
                message = "**** Exon nr " + str(microexon) + " was skipped to ease alignment"
                print(message)
                logging.info(message)
            # re-number exons (enumerate method starting from 1)
            #M = {}
            #for number, features in enumerate(Mm_exons, 1):
            #    M[number] = Mm_exons[features]
            #Mm_exons = M
        
        expected_exons = tuple(e for e, features in Mm_exons.items())
        message = "........................................................................\n "\
                "Expected exons : %s for protein : %s" % (str(expected_exons), protein_name)
        print(message)
        logging.info(message)
        
    return Mm_exons, expected_exons, microexons, microexons_seqs

"""