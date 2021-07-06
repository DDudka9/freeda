#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 18:23:12 2021

@author: damian

Uses the fact that dictionaries in Python 3 are ordered to iterate over an MSA
and finds exons, which are then identified as: missing, intronic, not-intronic or possibly
retrotransposition. It also calls synteny and duplications.

"""


# CONSIDER ELIMINATING FRAMESHIFR CHECK -> currently only detects insertions and not deletions
# use nomLeu3_genome (Ne) for CENPT -> exon 1 goes cose its 1bp deletion but exon 2 is missing cose 1bp insertion
# I ADDED DELETION COUNT TO INSERTION COUNT

import logging

def find_exons(cds, locus, gene, seqs, contig_name, Mm_exons, expected_exons, microexons): # works well

    exons = {}
    exon = ""
    exon_start = 0
    exon_end = 0
    exon_number = 0
    exon_checked = False
    introny = False
    introny_at_Nterm = False
    introny_at_Cterm = False
    exon_missing = False
    last_bp = False
    possible_retrotransposition = False
    synteny = False
    N_term_synteny = False
    C_term_synteny = False
    nr_of_intronic_exons = 0
    nr_of_RETRO_exons = 0
    RETRO_score = 0
    nr_of_fully_intronic_exons = 0
    nr_of_partially_intronic_exons = 0
    duplication_score = 0
    insertion = 0
    insertion_with_N = False
    big_insertion = False
    Cterm_synteny_message = None
    single_exon = True
    non_ACGT = False
    
    duplication_score_parameter = False # THIS WILL BE PART OF NEW MODULE
    
    # all sequences are in capital letters at this point
    list_of_non_ACGT = ["N","Y","R","W","S","K","M","D","H","V","B","X"]

    for position in cds:
        
        # NON-CODING REGION IN ORIGINAL SPECIES
                
        if exon_checked == False and cds[position] == "-":
            continue
        
        
        # EXON STARTS IN ORIGINAL SPECIES -> define N-term for contig locus
        
        if exon_checked == False and cds[position] == gene[position]:
            
            exon_number += 1
            
            # checkpoint for potentially skipped microexons
            if exon_number not in expected_exons:
                exon_number += 1

            exon = cds[position]
            exon_start = position
            exon_checked = True
            exon_missing = False
            N_term_retrotransposition = False
            C_term_retrotransposition = False
            divergent_introns = False
            
            # check if there is an exon in the analysed contig
            if locus[position] in "ACTG":

                # check if this exon is intronic at N-term
                if check_introny(position, last_bp, cds, locus, gene) == True:
                    introny_at_Nterm = True
                    
                    # check if its the first exon based on original species cds
                    if exon_number == 1:
                        
                        # check synteny
                        if check_synteny_Nterm(position, locus, gene) == True:
                            N_term_synteny = True

                        # do not allow introny in first exon if not syntenic (often RETRO have 5UTR)
                        else:
                            introny_at_Nterm = False
                
                # check possible retrotransposition since no intron attached
                elif locus[position-1] == "-" \
                    and locus[position-4] == "-" \
                    and locus[position-7] == "-":
                    
                    N_term_retrotransposition = check_retrotransposition(position, last_bp, cds, locus, gene)
            
                # check for very divergent introns (will be counted intronic)
                # does not allow the first exon to have divergent introns (possibly not-syntenic exon)
                elif exon_number != 1 and homology_check(position, last_bp, cds, locus, gene) == True:
                    divergent_introns = True
            
            # this exon seems to be missing from the contig          
            else:
                exon_missing = True
        
        # EXON CONTINUES IN ORIGINAL SPECIES (includes first bp)
    
        if exon_checked == True and cds[position] == gene[position]:
            
            # clone that position
            exon = exon + cds[position]
            
            # check for non_ACGT positions
            if locus[position] in list_of_non_ACGT:
                non_ACGT = True
            
            # count insertions in the analysed locus
            if cds[position] == "-" and gene[position] == "-" and locus[position] != "-":
                
                # Ns in insertions lead to misalignment
                if locus[position] == "N":
                    insertion_with_N = True
                    
                if insertion == 0:
                    insertion = 1
                    continue
                    #m = "insertion: %s %s" % (insertion, locus[position]) 
                    #print(m)
                    #logging.info(m)
                if insertion > 0:
                    insertion = insertion + 1
                    #m = "insertion: %s %s" % (insertion, locus[position])
                    #print(m)
                    #logging.info(m)
            
            # count deletions as well
            elif locus[position] == "-":
                
                if insertion == 0:
                    insertion = -1
                
                else:
                    insertion = insertion - 1
            
            
            # go to next position
            continue
        
        
        # EXON ENDS IN ORIGINAL SPECIES -> define C-term
        
        if exon_checked == True and cds[position] == "-" and gene[position] != "-":
            
            # mark start of the non-coding sequence
            last_bp = True
            exon_checked = False
            
            # check introny first
            if check_introny(position, last_bp, cds, locus, gene) == True:
                introny_at_Cterm = True
                    
            # check if its the last exon based on original species cds
            if exon_number == list(Mm_exons.keys())[-1]:
                        
                # check synteny
                # WARNING: average 800bp 3'UTRs in mammals make synteny check at C-term not very efficient (most contigs too short)
                C_term_synteny, Cterm_synteny_message = check_synteny_Cterm(position, locus, gene, contig_name, exon_number)
                if C_term_synteny == False:
                    
                    # do not allow introny in last exon if not syntenic (often RETRO have 5UTR)
                    introny_at_Cterm = False
            
            # check possible retrotransposition since no intron attached
            elif locus[position] == "-" \
                and locus[position+3] == "-" \
                and locus[position+6] == "-" and locus[position-1] not in "-N":
                                
                if check_retrotransposition(position, last_bp, cds, locus, gene) == True:
                    C_term_retrotransposition = True
            
            # check for very divergent introns (will be counted intronic) -> does not allow divergent introns in last exon
            elif introny_at_Cterm == False \
                and exon_number != list(Mm_exons.keys())[-1] \
                and homology_check(position, last_bp, cds, locus, gene) == True:
                divergent_introns = True
                
        # CALL EXON IN CONTIG LOCUS ANALYSED
        
        if last_bp == True:

            # call it missing
            if exon_missing == True:
                message = "Contig %s exon %s genomic locus is MISSING" \
                              % (contig_name, exon_number)
                print(message)
                logging.info(message)
                
                # this is the first bp of an intron -> will not be cloned  
                exon_end = position
                
                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, introny, exon_number, big_insertion
                
                # reset all parameters
                exon, insertion, introny_at_Nterm, introny_at_Cterm, \
                introny, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()
                
                # reset non_ACGT to False in case it was triggered
                non_ACGT = False
                
                continue


            # call it possible retrotransposition
            if N_term_retrotransposition == True or C_term_retrotransposition == True:
            
                # RETRO at both ends
                if N_term_retrotransposition == True and C_term_retrotransposition == True:
                    possible_retrotransposition = True
                    nr_of_RETRO_exons += 1
                    message = "Contig %s exon %s genomic locus might be RETRO (N and C-term)" \
                            % (contig_name, exon_number)
                    print(message)
                    logging.info(message)
            
                # RETRO at N-term
                elif N_term_retrotransposition == True and C_term_retrotransposition == False:
                    possible_retrotransposition = True
                    nr_of_RETRO_exons += 1
                    if introny_at_Cterm == True:
                        introny = True
                        
                        # check insertions
                        big_insertion = check_insertion(contig_name, insertion, exon_number, Mm_exons, insertion_with_N)
                        #message = "insertion: %s" % str(insertion)
                        #print(message)
                        #logging.info(message)
                        
                        message = "Contig %s exon %s genomic might be RETRO (N-term) but is intronic (C-term)" \
                            % (contig_name, exon_number)
                        print(message)
                        logging.info(message)
                    else:
                        message = "Contig %s exon %s genomic might be RETRO (N-term)" \
                            % (contig_name, exon_number)
                        print(message)
                        logging.info(message)
            
                # RETRO at C-term
                elif N_term_retrotransposition == False and C_term_retrotransposition == True:
                    possible_retrotransposition = True
                    nr_of_RETRO_exons += 1
                    if introny_at_Nterm == True:
                        introny = True
                        
                        # check insertions
                        big_insertion = check_insertion(contig_name, insertion, exon_number, Mm_exons, insertion_with_N)
                        #message = "insertion: %s" % str(insertion)
                        #print(message)
                        #logging.info(message)
                        
                        message = "Contig %s exon %s genomic is intronic (N-term) but might be RETRO (C-term)" \
                            % (contig_name, exon_number)
                        print(message)
                        logging.info(message)
                        
                    else:
                        message = "Contig %s exon %s genomic might be RETRO (C-term)" \
                            % (contig_name, exon_number)
                        print(message)
                        logging.info(message)
                
                # this is the first bp of an intron -> will not be cloned  
                exon_end = position
                
                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, introny, exon_number, big_insertion
                
                # reset all parameters
                exon, insertion, introny_at_Nterm, introny_at_Cterm, \
                introny, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()
                
                continue
            

            # call it non-intronic
            if (divergent_introns == False or (divergent_introns == True and single_exon == True)) \
                                    and introny_at_Nterm == False and introny_at_Cterm == False:
                message = "Contig %s exon %s genomic locus is not intronic" \
                        % (contig_name, exon_number)
                print(message)
                logging.info(message)
                
                # this is the first bp of an intron -> will not be cloned  
                exon_end = position
                
                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, introny, exon_number, big_insertion
                
                # reset all parameters
                exon, insertion, introny_at_Nterm, introny_at_Cterm, \
                introny, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()
                
                # reset non_ACGT to False in case it was triggered
                non_ACGT = False
                
                continue
            
            
            # call it invalid if non_ACGT present
            if non_ACGT == True:
                introny = False
                
                message = "Contig %s exon %s genomic locus contains non_ACGT bases (not allowed)" \
                        % (contig_name, exon_number)
                print(message)
                logging.info(message)
                
                # this is the first bp of an intron -> will not be cloned  
                exon_end = position
                
                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, introny, exon_number, big_insertion
                
                # reset all parameters
                exon, insertion, introny_at_Nterm, introny_at_Cterm, \
                introny, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()
                
                # reset non_ACGT parameter
                non_ACGT = False
                
                continue
            
            
            # call it divergent introns (only if its not a lone exon in contig -> risk of erronous introny)
            if divergent_introns == True and introny_at_Nterm == False and introny_at_Cterm == False \
                                                and single_exon == False:
                introny = True
                
                big_insertion = check_insertion(contig_name, insertion, exon_number, Mm_exons, insertion_with_N)
                #message = "insertion: %s" % str(insertion)
                #print(message)
                #logging.info(message)

                message = "Contig %s exon %s genomic locus has divergent introns (allowed)" \
                        % (contig_name, exon_number)
                print(message)
                logging.info(message)
                
                # add to pool of intronic exons
                nr_of_intronic_exons += 1
                
                # mark it as not being intronic at both ends
                nr_of_partially_intronic_exons += 1
                
                # this is the first bp of an intron -> will not be cloned  
                exon_end = position
                
                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, introny, exon_number, big_insertion
                
                # reset all parameters
                exon, insertion, introny_at_Nterm, introny_at_Cterm, \
                introny, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()
                
                continue
            
            
            # call it intronic
            if introny_at_Nterm == True or introny_at_Cterm == True:
                introny = True
                single_exon = False
                
                # at both ends
                if introny_at_Nterm == True and introny_at_Cterm == True:
                    message = "Contig %s exon %s genomic locus is INTRONIC  (N and C-term)" \
                        % (contig_name, exon_number)
                    print(message)
                    logging.info(message)   
                    
                    # mark it as fully intronic at both ends
                    nr_of_fully_intronic_exons += 1
                
                # at N-term
                elif introny_at_Nterm == True and introny_at_Cterm == False:
                    message = "Contig %s exon %s genomic locus is INTRONIC (only N-term)" \
                        % (contig_name, exon_number)
                    print(message)
                    logging.info(message)
                    
                    # mark it as not being intronic at both ends
                    nr_of_partially_intronic_exons += 1
                
                # at C-term
                elif introny_at_Nterm == False and introny_at_Cterm == True:
                    message = "Contig %s exon %s genomic locus is INTRONIC (only C-term)" \
                        % (contig_name, exon_number)
                    print(message)
                    logging.info(message)
                    
                    # mark it as not being intronic at both ends
                    nr_of_partially_intronic_exons += 1
                
                # CATCH NO CALL DECISIONS
                else:
                    message = "No call in introny: \nintrony = %s, \nN_term introny = %s, \nC_term introny = %s, \ndivergent exons = %s" \
                    "\nN_term_retrotransposition = %s, \nC_term_retrotransposition = %s" \
                    % (introny, introny_at_Nterm, introny_at_Cterm, divergent_introns, N_term_retrotransposition, C_term_retrotransposition)
                    print(message)
                    logging.info(message)
                    
                    
                # add to pool of intronic exons
                nr_of_intronic_exons += 1
                
                # check insertions
                big_insertion = check_insertion(contig_name, insertion, exon_number, Mm_exons, insertion_with_N)
                #message = "insertion: %s" % str(insertion)
                #print(message)
                #logging.info(message)
                
                # this is the first bp of an intron -> will not be cloned  
                exon_end = position
                
                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, introny, exon_number, big_insertion
                
                # reset all parameters
                exon, insertion, introny_at_Nterm, introny_at_Cterm, \
                introny, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()

                continue
            
            # CATCH NO CALL DECISIONS
            else:
                message = "No call at all: \nintrony = %s, \nN_term introny = %s, \nC_term introny = %s, \ndivergent introns = %s" \
                    "\nsingle exon = %s, \nN_term_retrotransposition = %s, \nC_term_retrotransposition = %s" \
                    % (introny, introny_at_Nterm, introny_at_Cterm, divergent_introns, single_exon, N_term_retrotransposition, C_term_retrotransposition)
                print(message)
                logging.info(message)
                
    
    print(Cterm_synteny_message)
    logging.info(Cterm_synteny_message)
    
    # check the probability of a false positive RETRO call for this contig
    if nr_of_RETRO_exons != 0 and nr_of_intronic_exons != 0:
        RETRO_score = float(nr_of_RETRO_exons/nr_of_intronic_exons)
    elif nr_of_RETRO_exons == 0:
        RETRO_score = float(0)
    elif nr_of_RETRO_exons != 0 and nr_of_intronic_exons == 0:
        RETRO_score = float(1)
    
    message = "RETRO_score (>= 0.4 suggests retrotransposition event) = %s" % (str(RETRO_score))
    print(message)
    logging.info(message)
    
    if RETRO_score < 0.4:
    
        
    # CALL SYNTENY OF CONTIG LOCUS ANALYSED
    
        # call synteny of the contig
        if N_term_synteny == True and C_term_synteny == True:
            message = "\nSYNTENY ESTIMATION -> Contig %s is SYNTENIC (N and C-term)" % (contig_name)
            print(message)
            logging.info(message)
            synteny = True
    
        elif N_term_synteny == True and C_term_synteny == False:
            message = "\nSYNTENY ESTIMATION -> Contig %s is SYNTENIC (N-term)" % (contig_name)
            print(message)
            logging.info(message)
            synteny = True
    
        elif N_term_synteny == False and C_term_synteny == True:
            message = "\nSYNTENY ESTIMATION -> Contig %s is SYNTENIC (C-term)" % (contig_name)
            print(message)
            logging.info(message)
            synteny = True
            
        elif N_term_synteny == False and C_term_synteny == False:
            message = "\nSYNTENY ESTIMATION -> Contig %s is NOT SYNTENIC" % (contig_name)
            print(message)
            logging.info(message)
    
    else:
        message = "\nSYNTENY ESTIMATION: possible retrotransposition (synteny not allowed)"
        print(message)
        logging.info(message)
        
    if exon_number != len(Mm_exons):
        message = "\nFREEDA could not find all exons in contig %s " % (contig_name)
        print(message)
        logging.info(message)
        
        if microexons != []:
            message = "\nFREEDA skipped microexons: %s " % (microexons)
            print(message)
            logging.info(message)
        
    # assess likelihood of this contig carrying a duplication
    if duplication_score_parameter == True:
        
        if nr_of_intronic_exons != 0 and nr_of_fully_intronic_exons != 0:
            duplication_score = nr_of_fully_intronic_exons/nr_of_intronic_exons
    
            message = "\nDuplication_score (<0.5 suggests duplication) = %s" % (str(duplication_score))
            print(message)
            logging.info(message)
            
    if duplication_score_parameter == False:
        duplication_score = 1.0
        
        message = "\nDuplication_score (<0.5 suggests duplication) = DISABLED (1.0 as default)"
        print(message)
        logging.info(message)

    
    return exons, possible_retrotransposition, synteny, RETRO_score, duplication_score


def check_introny(position, last_bp, cds, locus, gene): # works well
    
    # default state of introny is false
    introny = False
    
    max_length = len(locus)
    extension1 = 100
    extension2 = 50
    extension3 = 20
    
    min_stretch_length = 50
    homology_threshold = 0.75
    no_homology_threshold = 0.66
    allowed_indels = 20
    
    # for testing introny at the beginning of an exon
    if last_bp == False:
        # mark the original posisiton to be used for homology_check
        starting_position = position
        # start from previous position
        position -= 1
        offset = -1
    
        # define how far upstream can rolling hd go
        if position >= extension1 - 1:
            end = position - extension1
            extension = extension1
            
        elif position < extension1 - 1 and position >= extension2 - 1:
            end = position - extension2
            extension = extension2
            
        elif position < extension2 - 1 and position >= extension3 - 1:
            end = position - extension3
            extension = extension3
            
        else:
            return introny
    
    
    if last_bp == True:
        # mark the original posisiton to be used for homology_check
        starting_position = position
        offset = 1
        
        # define how far downstream can rolling hd go
        if position + extension1 <= max_length + 1:
            end = position + extension1 - 1
            extension = extension1
            
        elif position + extension1 > max_length + 1 and position + extension2 <= max_length + 1:
            end = position + extension2 - 1
            extension = extension2
            
        elif position + extension2 > max_length + 1 and position + extension3 <= max_length + 1:
            end = position + extension3 - 1
            extension = extension3
            
        else:
            return introny
    
    #if extension == 20 or extension == 50:
        #message = "\n$$$$ extension: %s\n" % extension
        #print(message)
        #logging.info(message)
    
    stretch_length = 0
    matches = 0
    indels = 0   
    
    
    for i in range(position, end, offset):
        
        if stretch_length != 0:
            rolling_hd = matches / stretch_length
            #message = "> dynamic hamming distance: %s, stretch: %s bp" % (rolling_hd, stretch_length)
            #print(message)
            #logging.info(message)
        
        # detect homology at minimum stretch 
        if stretch_length >= min_stretch_length and rolling_hd >= homology_threshold:
            #message = ">>> reached introny : %s, stretch: %s bp" % (rolling_hd, stretch_length)
            #print(message)
            #logging.info(message)
            introny = True
            return introny
        
        # detect lack of homology at minimum stretch 
        elif stretch_length >= min_stretch_length and rolling_hd <= no_homology_threshold:
            #message = ">>> failed introny : %s, stretch: %s bp" % (rolling_hd, stretch_length)
            #print(message)
            #logging.info(message)
            return introny
        
        # no homology until extension length
        elif stretch_length == extension and rolling_hd < homology_threshold:
            #message = ">>> failed introny : %s, stretch: %s bp" % (rolling_hd, stretch_length)
            #print(message)
            #logging.info(message)
            return introny
        
        # ADDED "cds[i] == "-":" at 5.30pm on 03/04/2021
        
        # count matches and mismatches along the extension
        elif stretch_length < extension and cds[i] == "-":
            
            if indels < allowed_indels:
                
                # count but dont yet add locus indels to stretch length
                if locus[i] == "-" or gene[i] == "-":
                    indels += 1
                    continue
                
                elif locus[i] == gene[i]:
                    matches += 1
                    stretch_length += 1
                
                else:
                    stretch_length += 1
                    
            if indels >= allowed_indels:
                
                # count and add locus indels to stretch length
                if locus[i] == "-" or gene[i] == "-":
                    stretch_length += 1
                    continue
                
                elif locus[i] == gene[i]:
                    matches += 1
                    stretch_length += 1

                else:
                    stretch_length += 1
            
        # if introny reached at the last bp of the extension
        elif stretch_length == extension and rolling_hd >= homology_threshold:
            
            # check homology inwards before calling introny
            homology = homology_check(starting_position, last_bp, cds, locus, gene)
            if homology == True:
                introny = True
                return introny
            else:
                return introny
        
        # the introny check function clashed with another exon -> no introny as default
        # sometimes single bp would be coopted from introns and be treated as start of another exon
        elif cds[i] != "-" and cds[i-1] != "-":
            message = "introny = %s ------ another exon interferes with introny check -> no introny" % introny
            print(message)
            logging.info(message)
            return introny
        
        # last bp and no introny
        elif i == end-1 and rolling_hd < homology_threshold:
            return introny
        
        # last bp and introny
        elif i == end-1 and rolling_hd >= homology_threshold:
            introny = True
            return introny
            
    
    # loop ended and introny wasnt called -> default False 
    return introny
    

def homology_check(starting_position, last_bp, cds, locus, gene): # runs only when ambigous introny (100bp stretch with 0.80 homology) or RETRO suspected
    
    homology = False
    
    min_stretch_length = 30
    max_stretch_length = 100
    homology_threshold = 0.80
    no_homology_threshold = 0.66
    stretch_length = 0
    matches = 0
    indels = 0
    allowed_indels = 10
    
    if last_bp == False:
        
        # cannot predict the exon length at this point so hardcoded 100 as stop
        for i in range(starting_position, starting_position + max_stretch_length):
            
            if stretch_length < min_stretch_length and indels < allowed_indels:
                
                #message = "* N - stretch : %s bp and matches : %s" % (stretch_length, matches)
                #logging.info(message)
                #print(message)
            
                if locus[i] == cds[i]:
                    matches += 1
                    stretch_length += 1
                    continue
                
                # dont count indels within min stretch
                elif locus[i] == "-" or cds[i] == "-":
                    indels += 1
                    continue
                
                # rare case of homology check going beyond the exon in question
                elif locus[i] == "-" or gene[i] == "-":
                    indels += 1
                    continue
                
                # mismatch
                else:
                    stretch_length += 1
                    continue
            
            if stretch_length >= min_stretch_length:
                
                # estimate homology
                rolling_hd = matches / stretch_length
                #message = "*** N - homology : %s and stretch: %s bp" % (rolling_hd, stretch_length)
                #print(message)
                #logging.info(message)
            
                # call homology if passed min_stretch_length
                if rolling_hd >= homology_threshold:
                    homology = True
                    #message = "****** N - HOMOLOGY : %s and stretch: %s bp" % (rolling_hd, stretch_length)
                    #print(message)
                    #logging.info(message)
                    return homology
        
                # call lack of homology if passed min_stretch_length and no_homology_threshold is reached
                elif rolling_hd <= no_homology_threshold:
                    #message = "****** N - NO homology : %s and stretch: %s bp" % (rolling_hd, stretch_length)
                    #print(message)
                    #logging.info(message)
                    return homology
        
                # call lack of homology if reached max_stretch_length but not the homology threshold
                elif stretch_length == max_stretch_length and rolling_hd < homology_threshold:
                    #message = "****** N - NO homology : %s and stretch: %s bp" % (rolling_hd, stretch_length)
                    #print(message)
                    #logging.info(message)
                    return homology
                
                # not ready to call homology
                elif rolling_hd < homology_threshold and stretch_length < max_stretch_length:
                    if locus[i] == cds[i]:
                        matches += 1
                        stretch_length += 1
                    
                    # mismatch or indel
                    if locus[i] != cds[i]:
                        stretch_length += 1
            
            # too many indels with too little bp aligned (RISK OF LOOSING RETRO EXONS WITH PERIPHERAL INDELS)
            if indels >= allowed_indels and stretch_length < min_stretch_length:
                return homology
        
            # end of the exon
            if cds[i] != gene[i]:
                return homology
        
        # if loop is finished and homology not called -> default False 
        return homology
        

    if last_bp == True:
        # set position to last exon bp
        starting_position -= 1
        
        # cannot predict the exon length at this point so hardcoded 100 as stop
        for i in range(starting_position, starting_position - max_stretch_length, -1):
            
            if stretch_length < min_stretch_length and indels < allowed_indels:
                
                #message = "* C - stretch : %s bp and matches : %s" % (stretch_length, matches)
                #print(message)
                #logging.info(message)
            
                if locus[i] == cds[i]:
                    matches += 1
                    stretch_length += 1
                    continue
                
                # dont count indels within min stretch
                elif locus[i] == "-" or cds[i] == "-":
                    indels += 1
                    continue
                
                # rare case of homology check going beyond the exon in question
                elif locus[i] == "-" or gene[i] == "-":
                    indels += 1
                    continue
                
                # mismatch
                else:
                    stretch_length += 1
                    continue
            
            if stretch_length >= min_stretch_length:
                
                # estimate homology
                rolling_hd = matches / stretch_length
                #message = "*** C - homology : %s and stretch: %s bp" % (rolling_hd, stretch_length)
                #print(message)
                #logging.info(message)
            
                # call homology if passed min_stretch_length
                if rolling_hd >= homology_threshold:
                    #message = "****** C - HOMOLOGY : %s and stretch: %s bp" % (rolling_hd, stretch_length)
                    #print(message)
                    #logging.info(message)
                    homology = True
                    return homology
        
                # call lack of homology if passed min_stretch_length and no_homology_threshold is reached
                elif rolling_hd <= no_homology_threshold:
                    #message = "****** C - NO homology : %s and stretch: %s bp" % (rolling_hd, stretch_length)
                    #print(message)
                    #logging.info(message)
                    return homology
        
                # call lack of homology if reached max_stretch_length but not the homology threshold
                elif stretch_length == max_stretch_length and rolling_hd < homology_threshold:
                    #message = "****** C- NO homology : %s and stretch: %s bp" % (rolling_hd, stretch_length)
                    #print(message)
                    #logging.info(message)
                    return homology
                
                # not ready to call homology
                elif stretch_length < max_stretch_length and rolling_hd < homology_threshold:
                    if locus[i] == cds[i]:
                        matches += 1
                        stretch_length += 1
                    
                    # mismatch or indel
                    if locus[i] != cds[i]:
                        stretch_length += 1
            
            # too many indels with too little bp aligned (RISK OF LOOSING RETRO EXONS WITH PERIPHERAL INDELS)
            if indels >= allowed_indels and stretch_length < min_stretch_length:
                return homology
            
            # end of the exon
            if cds[i] != gene[i]:
                return homology
            
        # if loop is finished and homology not called -> default False        
        return homology

    
def check_synteny_Nterm(starting_position, locus, gene):
    
    synteny = False
    UTR_length = 200
    stretch_length = 0
    matches = 0
    homology_threshold = 0.75
    
    # start from previous position
    starting_position -= 1
    
    # if contig too short, call lack of synteny
    if starting_position - UTR_length < 0:
        
        message = "\n..............Contig too short to check N-term synteny \n"
        print(message)
        logging.info(message)
        
        return synteny
    
    # if contig is long enough, estimate synteny
    for i in range(starting_position, starting_position - UTR_length + 1, -1):
        
        # skip indels
        if locus[i] == "-" or gene[i] == "-":
            continue
        
        # match
        elif locus[i] == gene[i]:
            matches += 1
            stretch_length += 1
        
        # mismatch
        elif locus[i] != gene[i]:
            stretch_length += 1
    
    # avoid divisions over 0
    if stretch_length == 0:
        
        message = "()()() stretch_length == 0"
        print(message)
        logging.info(message)
        
        return synteny
    
    # cannot call synteny if too many indels
    if stretch_length < 100:
        
        message = "\n..............Contig too short to check N-term synteny: only %s bp aligned \n" % str(stretch_length)
        print(message)
        logging.info(message)
        
        return synteny
    
    # check synteny
    if matches/stretch_length >= homology_threshold:
        synteny = True
        message = "\n..............N-term syntenic: %s hamming distance\n" % str(matches/stretch_length)
        print(message)
        logging.info(message)
        
    else:
        synteny = False
        message = "\n..............N-term NOT syntenic: %s hamming distance\n" % str(matches/stretch_length)
        print(message)
        logging.info(message)
    
    return synteny


def check_synteny_Cterm(starting_position, locus, gene, contig_name, exon_number):
    
    synteny = False
    UTR_length = 200
    stretch_length = 0
    matches = 0
    homology_threshold = 0.75
    
    # if contig too short, call lack of synteny
    if starting_position + UTR_length > len(locus):
        
        Cterm_synteny_message = "\n.............Alignment is too short to check C-term synteny \n"        
        return synteny, Cterm_synteny_message
    
    # if contig is long enough, estimate synteny
    for i in range(starting_position, starting_position + UTR_length):
        
        # skip indels
        if locus[i] == "-" or gene[i] == "-":
            continue
        
        # match
        elif locus[i] == gene[i]:
            matches += 1
            stretch_length += 1
        
        # mismatch
        elif locus[i] != gene[i]:
            stretch_length += 1

    
    # avoid divisions over 0
    if stretch_length == 0:
        
        Cterm_synteny_message = "\n.............Contig too short to check C-term synteny: only %s bp aligned\n" \
                                                           % str(stretch_length)      
        return synteny, Cterm_synteny_message
    
    # cannot call synteny if too many indels
    if stretch_length < 100:
        
        Cterm_synteny_message = "\n.............Contig too short to check C-term synteny: only %s bp aligned\n" % str(stretch_length)
        return synteny, Cterm_synteny_message
        
    # check synteny
    if matches/stretch_length >= homology_threshold:
        synteny = True
        Cterm_synteny_message = "\n..............C-term syntenic: %s hamming distance\n" % str(matches/stretch_length)
        
    else:
        synteny = False
        Cterm_synteny_message = "\n..............C-term NOT syntenic: %s hamming distance\n" % str(matches/stretch_length)
    
    return synteny, Cterm_synteny_message


def check_retrotransposition(position, last_bp, cds, locus, gene): # runs only when introny is False
    
    retrotransposition = False
    
    if homology_check(position, last_bp, cds, locus, gene) == True:
        retrotransposition = True
        return retrotransposition
    else:
        return retrotransposition


def reset_exon_parameters():
    
     exon = ""
     insertion = 0
     introny_at_Nterm = False
     introny_at_Cterm = False
     introny = False
     exon_missing = False
     last_bp = False
     big_insertion = False
     insertion_with_N = False
     
     return exon, insertion, introny_at_Nterm, introny_at_Cterm, \
         introny, exon_missing, last_bp, big_insertion, insertion_with_N


def check_insertion(contig_name, insertion, exon_number, Mm_exons, insertion_with_N):
    big_insertion = False
    
    # do not allow Ns in insertions to avoid misalignment
    if insertion_with_N == True:
        big_insertion = True
        insertion_message = ">>>>> Contig %s exon %s contains %s bp insertion (Ns present -> not allowed) <<<<< \n" \
            % (contig_name, exon_number, str(insertion))
        print(insertion_message)
        logging.info(insertion_message)
        return big_insertion
    
    # allow insertions (even frameshifts) in last exons
    if exon_number == list(Mm_exons.keys())[-1] and (insertion >= 180 or (insertion >= 10 and insertion % 3 != 0)):
        insertion_message = ">>>>> Contig %s exon %s contains %s bp insertion (LAST EXON -> allowed) <<<<< \n" \
            % (contig_name, exon_number, str(insertion))
        print(insertion_message)
        logging.info(insertion_message)
        return big_insertion
    
    # do not allow big insertions or frameshift insertions in middle exons (conservative approach)
    if insertion >= 180: # or (insertion >= 10 and insertion % 3 != 0):
        big_insertion = True
        insertion_message = ">>>>> Contig %s exon %s contains %s bp insertion (not allowed) <<<<< \n" \
            % (contig_name, exon_number, str(insertion))
        print(insertion_message)
        logging.info(insertion_message)
        return big_insertion
    
    else:
        return big_insertion
