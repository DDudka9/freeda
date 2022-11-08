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

Finds exons in genomic assemblies

"""

import logging


def find_exons(gene_name, cds_seq, locus_seq, gene_seq, contig_name, ref_exons, expected_exons, all_genes_dict=None):
    """Finds and calls exons based on the cds_seq and exon make-up from reference species"""

    last_exon = expected_exons[-1]
    exons = {}
    exon = ""
    exon_start = 0
    exon_number = min(list(map(int, expected_exons))) - 1  # gets the smallest exon number expected - 1 (usually 0)
    exon_checked = False
    intron = False
    intron_at_5_prime = False
    intron_at_3_prime = False
    exon_missing = False
    last_bp = False
    possible_retrotransposition = False
    synteny = False
    five_prime_synteny = False
    three_prime_synteny = False
    nr_of_intron_exons = 0
    nr_of_RETRO_exons = 0
    RETRO_score = 0
    nr_of_fully_intron_exons = 0
    nr_of_partially_intron_exons = 0
    synteny_score = 0
    insertion = 0
    insertion_with_N = False
    big_insertion = False
    three_prime_synteny_message = None
    single_exon = True
    non_ACGT = False
    synteny_score_parameter = False

    if all_genes_dict:
        if all_genes_dict[gene_name][0] == "Duplication expected":
            synteny_score_parameter = True

    # all sequences are in capital letters at this point
    list_of_non_ACGT = ["N", "Y", "R", "W", "S", "K", "M", "D", "H", "V", "B", "X"]

    for position in cds_seq:


        # NON-CODING REGION IN REF SPECIES

        if exon_checked is False and cds_seq[position] == "-":
            continue


        # EXON STARTS IN REF SPECIES -> define 5-prime for contig locus_seq

        if exon_checked is False and cds_seq[position] == gene_seq[position]:

            exon_number += 1

            # checkpoint for potentially skipped microexons
            if exon_number not in expected_exons:
                exon_number += 1

            exon = cds_seq[position]
            exon_start = position
            exon_checked = True
            exon_missing = False
            five_prime_retrotransposition = False
            three_prime_retrotransposition = False
            divergent_introns = False

            # check if there is an exon in the analyzed contig
            if locus_seq[position] in "ACTG":

                # check if this exon has intron at 5-prime
                if check_intron(position, last_bp, cds_seq, locus_seq, gene_seq) is True:
                    intron_at_5_prime = True

                    # check if its the first exon based on ref species cds_seq
                    if exon_number == 1:

                        # check synteny
                        if check_synteny_5_prime(position, locus_seq, gene_seq) is True:
                            five_prime_synteny = True

                        # do not allow intron in first exon if not syntenic (often RETRO have 5UTR)
                        else:
                            intron_at_5_prime = False

                else:
                    try:
                        # check possible retrotransposition since no intron attached
                        if locus_seq[position-1] == "-" \
                                and locus_seq[position-4] == "-" \
                                and locus_seq[position-7] == "-":

                            five_prime_retrotransposition = check_retrotransposition(position, last_bp, cds_seq,
                                                                                     locus_seq, gene_seq)

                    # KeyError triggered when match at the edge of alignment)
                    except KeyError:
                        intron_at_5_prime = False

            # this exon seems to be missing from the contig
            else:
                exon_missing = True


        # EXON CONTINUES IN REF SPECIES (includes first bp)

        if exon_checked is True and cds_seq[position] == gene_seq[position]:

            # clone that position
            exon = exon + cds_seq[position]

            # check for non_ACGT positions
            if locus_seq[position] in list_of_non_ACGT:
                non_ACGT = True

            # count insertions in the analyzed locus_seq
            if cds_seq[position] == "-" and gene_seq[position] == "-" and locus_seq[position] != "-":

                # Ns in insertions lead to misalignment
                if locus_seq[position] == "N":
                    insertion_with_N = True

                if insertion == 0:
                    insertion = 1
                    continue

                if insertion > 0:
                    insertion = insertion + 1

            # count deletions as well
            elif locus_seq[position] == "-":

                if insertion == 0:
                    insertion = -1

                else:
                    insertion = insertion - 1

            # go to next position
            continue


        # EXON ENDS IN REF SPECIES -> define 3_prime

        if exon_checked is True and cds_seq[position] == "-" and gene_seq[position] != "-":

            # mark start of the non-coding sequence
            last_bp = True
            exon_checked = False

            # check intron first
            if check_intron(position, last_bp, cds_seq, locus_seq, gene_seq) is True:
                intron_at_3_prime = True

            # check if its the last exon based on ref species cds_seq
            if exon_number == last_exon:

                # check synteny
                # WARNING: average 800bp 3'UTRs in mammals make synteny check at 3_prime not very efficient
                # (many contigs too short)
                three_prime_synteny, three_prime_synteny_message = check_synteny_3_prime(position, locus_seq, gene_seq)
                if three_prime_synteny is False:

                    # do not allow intron in last exon if not syntenic (often RETRO have 5UTR)
                    intron_at_3_prime = False

            # check possible retrotransposition since no intron attached -> retro is never checked in the last exon
            # but last exon retro seem to always cary a chunk of 3'UTR so last exon would unlikely ever be retro
            else:
                try:
                    if locus_seq[position] == "-" \
                                and locus_seq[position+3] == "-" \
                                and locus_seq[position+6] == "-" \
                                and locus_seq[position-1] not in "-N":

                        if check_retrotransposition(position, last_bp, cds_seq, locus_seq, gene_seq) is True:
                            three_prime_retrotransposition = True
                # KeyError triggered when match at the edge of alignment)
                except KeyError:
                    intron_at_3_prime = False

                # check for very divergent introns (will be counted as intron)
                # does not check for divergent introns in last exon
                # it will not run if exception is raised
                else:
                    if intron_at_3_prime is False \
                        and three_prime_retrotransposition is False \
                        and homology_check(position, last_bp, cds_seq, locus_seq, gene_seq) is True:
                        divergent_introns = True


        # CALL EXON IN CONTIG LOCUS ANALYZED

        if last_bp is True:

            # call it missing
            if exon_missing is True:
                message = "Contig %s exon %s genomic locus is MISSING" \
                              % (contig_name, exon_number)
                print(message)
                logging.info(message)

                # this is the first bp of an intron -> will not be cloned  
                exon_end = position

                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, intron, exon_number, big_insertion

                # reset all parameters
                exon, insertion, intron_at_5_prime, intron_at_3_prime, \
                intron, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()

                # reset non_ACGT to False in case it was triggered
                non_ACGT = False

                continue

            # call it possible retrotransposition
            if five_prime_retrotransposition is True or three_prime_retrotransposition is True:

                # RETRO at both ends
                if five_prime_retrotransposition is True and three_prime_retrotransposition is True:
                    possible_retrotransposition = True
                    nr_of_RETRO_exons += 1
                    message = "Contig %s exon %s genomic locus might be RETRO (5 and 3-prime)" \
                            % (contig_name, exon_number)
                    print(message)
                    logging.info(message)

                # RETRO at 5-prime
                elif five_prime_retrotransposition is True and three_prime_retrotransposition is False:
                    possible_retrotransposition = True
                    nr_of_RETRO_exons += 1
                    if intron_at_3_prime is True:
                        intron = True

                        # check insertions
                        big_insertion = check_insertion(contig_name, insertion, exon_number, ref_exons, insertion_with_N)

                        message = "Contig %s exon %s genomic might be RETRO (5-prime) but has intron (3-prime)" \
                            % (contig_name, exon_number)
                        print(message)
                        logging.info(message)
                    else:
                        message = "Contig %s exon %s genomic might be RETRO (5-prime)" \
                            % (contig_name, exon_number)
                        print(message)
                        logging.info(message)

                # RETRO at 3-prime
                elif five_prime_retrotransposition is False and three_prime_retrotransposition is True:
                    possible_retrotransposition = True
                    nr_of_RETRO_exons += 1
                    if intron_at_5_prime is True:
                        intron = True

                        # check insertions
                        big_insertion = check_insertion(contig_name, insertion, exon_number, ref_exons, insertion_with_N)

                        message = "Contig %s exon %s genomic has intron (5-prime) but might be RETRO (3-prime)" \
                            % (contig_name, exon_number)
                        print(message)
                        logging.info(message)

                    else:
                        message = "Contig %s exon %s genomic might be RETRO (3-prime)" \
                            % (contig_name, exon_number)
                        print(message)
                        logging.info(message)

                # this is the first bp of an intron -> will not be cloned  
                exon_end = position

                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, intron, exon_number, big_insertion

                # reset all parameters
                exon, insertion, intron_at_5_prime, intron_at_3_prime, \
                intron, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()

                continue

            # call it non-intron
            if (divergent_introns is False or (divergent_introns is True and single_exon is True)) \
                                    and intron_at_5_prime is False and intron_at_3_prime is False:
                message = "Contig %s exon %s genomic locus does not have intron" \
                        % (contig_name, exon_number)
                print(message)
                logging.info(message)

                # this is the first bp of an intron -> will not be cloned  
                exon_end = position

                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, intron, exon_number, big_insertion

                # reset all parameters
                exon, insertion, intron_at_5_prime, intron_at_3_prime, \
                intron, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()

                # reset non_ACGT to False in case it was triggered
                non_ACGT = False

                continue

            # call it invalid if non_ACGT present
            if non_ACGT is True:
                intron = False

                message = "Contig %s exon %s genomic locus contains non_ACGT bases (not allowed)" \
                        % (contig_name, exon_number)
                print(message)
                logging.info(message)

                # this is the first bp of an intron -> will not be cloned  
                exon_end = position

                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, intron, exon_number, big_insertion

                # reset all parameters
                exon, insertion, intron_at_5_prime, intron_at_3_prime, \
                intron, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()

                # reset non_ACGT parameter
                non_ACGT = False

                continue

            # call it divergent introns (only if its not a lone exon in contig -> risk of erronous intron)
            if divergent_introns is True and intron_at_5_prime is False and intron_at_3_prime is False \
                                                and single_exon is False:
                intron = True

                big_insertion = check_insertion(contig_name, insertion, exon_number, ref_exons, insertion_with_N)

                message = "Contig %s exon %s genomic locus has divergent introns (allowed)" \
                        % (contig_name, exon_number)
                print(message)
                logging.info(message)

                # add to pool of intron exons
                nr_of_intron_exons += 1

                # mark it as not being intron at both ends
                nr_of_partially_intron_exons += 1

                # this is the first bp of an intron -> will not be cloned  
                exon_end = position

                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, intron, exon_number, big_insertion

                # reset all parameters
                exon, insertion, intron_at_5_prime, intron_at_3_prime, \
                intron, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()

                continue

            # call it intron
            if intron_at_5_prime is True or intron_at_3_prime is True:
                intron = True
                single_exon = False

                # at both ends
                if intron_at_5_prime is True and intron_at_3_prime is True:
                    message = "Contig %s exon %s genomic locus has INTRON (5 and 3-prime)" \
                        % (contig_name, exon_number)
                    print(message)
                    logging.info(message)

                    # mark it as fully intron at both ends
                    nr_of_fully_intron_exons += 1

                # at 5-prime
                elif intron_at_5_prime is True and intron_at_3_prime is False:
                    message = "Contig %s exon %s genomic locus has INTRON (only 5-prime)" \
                        % (contig_name, exon_number)
                    print(message)
                    logging.info(message)

                    # mark it as not being intron at both ends
                    nr_of_partially_intron_exons += 1

                # at 3-prime
                elif intron_at_5_prime is False and intron_at_3_prime is True:
                    message = "Contig %s exon %s genomic locus has INTRON (only 3-prime)" \
                        % (contig_name, exon_number)
                    print(message)
                    logging.info(message)

                    # mark it as not being intron at both ends
                    nr_of_partially_intron_exons += 1

                # CATCH NO CALL DECISIONS
                else:
                    message = "No call in intron: \nintron = %s, \n5_prime intron = %s, \n3_prime intron = %s, " \
                              "\ndivergent exons = %s \n5_prime_retrotransposition = %s, \n3_prime_retrotransposition" \
                              " = %s" % (intron, intron_at_5_prime, intron_at_3_prime, divergent_introns,
                                         five_prime_retrotransposition, three_prime_retrotransposition)
                    print(message)
                    logging.info(message)


                # add to pool of intron exons
                nr_of_intron_exons += 1

                # check insertions
                big_insertion = check_insertion(contig_name, insertion, exon_number, ref_exons, insertion_with_N)

                # this is the first bp of an intron -> will not be cloned  
                exon_end = position

                # populate the exons dictionary with exon and its start:end as tuple
                exons[exon] = exon_start, exon_end, intron, exon_number, big_insertion

                # reset all parameters
                exon, insertion, intron_at_5_prime, intron_at_3_prime, \
                intron, exon_missing, last_bp, big_insertion, insertion_with_N = reset_exon_parameters()

                continue


            # CATCH NO CALL DECISIONS
            else:
                message = "No call at all: \nintron = %s, \n5_prime intron = %s, \n3_prime intron = %s, " \
                          "\ndivergent introns = %s" \
                    "\nsingle exon = %s, \n5_prime_retrotransposition = %s, \n3_prime_retrotransposition = %s" \
                    % (intron, intron_at_5_prime, intron_at_3_prime, divergent_introns, single_exon,
                       five_prime_retrotransposition, three_prime_retrotransposition)
                print(message)
                logging.info(message)

    print(three_prime_synteny_message)
    logging.info(three_prime_synteny_message)

    # check the probability of a false positive RETRO call for this contig
    if nr_of_RETRO_exons != 0 and nr_of_intron_exons != 0:
        RETRO_score = float(nr_of_RETRO_exons/nr_of_intron_exons)
    elif nr_of_RETRO_exons == 0:
        RETRO_score = float(0)
    elif nr_of_RETRO_exons != 0 and nr_of_intron_exons == 0:
        RETRO_score = float(100)

    message = "       RETRO_score (>= 0.4 suggests retrotransposition event) = %s" % (str(RETRO_score))
    print(message)
    logging.info(message)

    if RETRO_score < 0.4:


    # CALL SYNTENY OF CONTIG LOCUS ANALYZED

        # call synteny of the contig
        if five_prime_synteny is True and three_prime_synteny is True:
            message = "\n       SYNTENY ESTIMATION -> Contig %s is SYNTENIC (5 and 3-prime)" % (contig_name)
            print(message)
            logging.info(message)
            synteny = True

        elif five_prime_synteny is True and three_prime_synteny is False:
            message = "\n       SYNTENY ESTIMATION -> Contig %s is SYNTENIC (5-prime)" % (contig_name)
            print(message)
            logging.info(message)
            synteny = True

        elif five_prime_synteny is False and three_prime_synteny is True:
            message = "\n       SYNTENY ESTIMATION -> Contig %s is SYNTENIC (3-prime)" % (contig_name)
            print(message)
            logging.info(message)
            synteny = True

        elif five_prime_synteny is False and three_prime_synteny is False:
            message = "\n       SYNTENY ESTIMATION -> Contig %s is NOT SYNTENIC" % (contig_name)
            print(message)
            logging.info(message)

    else:
        message = "\n       SYNTENY ESTIMATION: possible retrotransposition (synteny not allowed)"
        print(message)
        logging.info(message)

    if exon_number != len(ref_exons): # removed "+ len(microexons)" 08_26_2021 (cose microexons are list of tuples)
        message = "\n       FREEDA could not find all exons in contig %s " % (contig_name)
        print(message)
        logging.info(message)

    # assess likelihood of this contig to be syntenic
    if synteny_score_parameter is True:

        if nr_of_intron_exons != 0 and nr_of_fully_intron_exons != 0:
            synteny_score = nr_of_fully_intron_exons/nr_of_intron_exons

            message = "\n       Synteny_score (<0.5 suggests poor synteny) = %s" % (str(synteny_score))
            print(message)
            logging.info(message)

    if synteny_score_parameter is False:
        synteny_score = 1.0

        message = "\n       Synteny_score (<0.5 suggests poor synteny) = DISABLED (1.0 as default)"
        print(message)
        logging.info(message)

    return exons, possible_retrotransposition, synteny, RETRO_score, synteny_score


def check_intron(position, last_bp, cds_seq, locus_seq, gene_seq):
    """Assesses if intron is present based on hamming distance scores"""

    # default state of intron is false
    intron = False

    max_length = len(locus_seq)
    extension1 = 100
    extension2 = 50
    extension3 = 20

    min_stretch_length = 50
    homology_threshold = 0.75
    no_homology_threshold = 0.66
    allowed_indels = 20

    # for testing intron at the beginning of an exon
    if last_bp is False:
        # mark the ref position to be used for homology_check
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
            return intron


    if last_bp is True:
        # mark the ref position to be used for homology_check
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
            return intron

    stretch_length = 0
    matches = 0
    indels = 0

    for i in range(position, end, offset):

        if stretch_length != 0:
            rolling_hd = matches / stretch_length

        # detect homology at minimum stretch 
        if stretch_length >= min_stretch_length and rolling_hd >= homology_threshold:
            intron = True
            return intron

        # detect lack of homology at minimum stretch 
        elif stretch_length >= min_stretch_length and rolling_hd <= no_homology_threshold:
            return intron

        # no homology until extension length
        elif stretch_length == extension and rolling_hd < homology_threshold:
            return intron

        # count matches and mismatches along the extension
        elif stretch_length < extension and cds_seq[i] == "-":

            if indels < allowed_indels:

                # count but dont yet add locus_seq indels to stretch length
                if locus_seq[i] == "-" or gene_seq[i] == "-":
                    indels += 1
                    continue

                elif locus_seq[i] == gene_seq[i]:
                    matches += 1
                    stretch_length += 1

                else:
                    stretch_length += 1

            if indels >= allowed_indels:

                # count and add locus_seq indels to stretch length
                if locus_seq[i] == "-" or gene_seq[i] == "-":
                    stretch_length += 1
                    continue

                elif locus_seq[i] == gene_seq[i]:
                    matches += 1
                    stretch_length += 1

                else:
                    stretch_length += 1

        # if intron reached at the last bp of the extension
        elif stretch_length == extension and rolling_hd >= homology_threshold:

            # check homology inwards before calling intron
            homology = homology_check(starting_position, last_bp, cds_seq, locus_seq, gene_seq)
            if homology is True:
                intron = True
                return intron
            else:
                return intron

        # the intron check function clashed with another exon -> no intron as default
        # sometimes single bp would be coopted from introns and be treated as start of another exon
        elif cds_seq[i] != "-" and cds_seq[i-1] != "-":
            message = "intron = %s ------ another exon interferes with intron check -> no intron" % intron
            print(message)
            logging.info(message)
            return intron

        # last bp and no intron
        elif i == end-1 and rolling_hd < homology_threshold:
            return intron

        # last bp and intron
        elif i == end-1 and rolling_hd >= homology_threshold:
            intron = True
            return intron

    # loop ended and intron wasnt called -> default False
    return intron


def homology_check(starting_position, last_bp, cds_seq, locus_seq, gene_seq):
    """Checks exon homology if: ambigous intron, RETRO suspected or truncation suspected at 5-prime"""
    # runs only when ambigous intron (100bp stretch with 0.80 homology) or RETRO suspected

    homology = False

    min_stretch_length = 30
    max_stretch_length = 100
    homology_threshold = 0.80
    no_homology_threshold = 0.66
    stretch_length = 0
    matches = 0
    indels = 0
    allowed_indels = 10

    if last_bp is False:

        # cannot predict the exon length at this point so hardcoded 100 as stop
        for i in range(starting_position, starting_position + max_stretch_length):

            if stretch_length < min_stretch_length and indels < allowed_indels:

                if locus_seq[i] == cds_seq[i]:
                    matches += 1
                    stretch_length += 1
                    continue

                # dont count indels within min stretch
                elif locus_seq[i] == "-" or cds_seq[i] == "-":
                    indels += 1
                    continue

                # rare case of homology check going beyond the exon in question
                elif locus_seq[i] == "-" or gene_seq[i] == "-":
                    indels += 1
                    continue

                # mismatch
                else:
                    stretch_length += 1
                    continue

            if stretch_length >= min_stretch_length:

                # estimate homology
                rolling_hd = matches / stretch_length

                # call homology if passed min_stretch_length
                if rolling_hd >= homology_threshold:
                    homology = True
                    return homology

                # call lack of homology if passed min_stretch_length and no_homology_threshold is reached
                elif rolling_hd <= no_homology_threshold:
                    return homology

                # call lack of homology if reached max_stretch_length but not the homology threshold
                elif stretch_length == max_stretch_length and rolling_hd < homology_threshold:
                    return homology

                # not ready to call homology
                elif rolling_hd < homology_threshold and stretch_length < max_stretch_length:
                    if locus_seq[i] == cds_seq[i]:
                        matches += 1
                        stretch_length += 1

                    # mismatch or indel
                    if locus_seq[i] != cds_seq[i]:
                        stretch_length += 1

            # too many indels with too little bp aligned (RISK OF LOOSING RETRO EXONS WITH PERIPHERAL INDELS)
            if indels >= allowed_indels and stretch_length < min_stretch_length:
                return homology

            # end of the exon
            if cds_seq[i] != gene_seq[i]:
                return homology

        # if loop is finished and homology not called -> default False 
        return homology


    if last_bp is True:
        # set position to last exon bp
        starting_position -= 1

        # cannot predict the exon length at this point so hardcoded 100 as stop
        for i in range(starting_position, starting_position - max_stretch_length, -1):

            if stretch_length < min_stretch_length and indels < allowed_indels:

                if locus_seq[i] == cds_seq[i]:
                    matches += 1
                    stretch_length += 1
                    continue

                # dont count indels within min stretch
                elif locus_seq[i] == "-" or cds_seq[i] == "-":
                    indels += 1
                    continue

                # rare case of homology check going beyond the exon in question
                elif locus_seq[i] == "-" or gene_seq[i] == "-":
                    indels += 1
                    continue

                # mismatch
                else:
                    stretch_length += 1
                    continue

            if stretch_length >= min_stretch_length:

                # estimate homology
                rolling_hd = matches / stretch_length

                # call homology if passed min_stretch_length
                if rolling_hd >= homology_threshold:
                    homology = True
                    return homology

                # call lack of homology if passed min_stretch_length and no_homology_threshold is reached
                elif rolling_hd <= no_homology_threshold:
                    return homology

                # call lack of homology if reached max_stretch_length but not the homology threshold
                elif stretch_length == max_stretch_length and rolling_hd < homology_threshold:
                    return homology

                # not ready to call homology
                elif stretch_length < max_stretch_length and rolling_hd < homology_threshold:
                    if locus_seq[i] == cds_seq[i]:
                        matches += 1
                        stretch_length += 1

                    # mismatch or indel
                    if locus_seq[i] != cds_seq[i]:
                        stretch_length += 1

            # too many indels with too little bp aligned (RISK OF LOOSING RETRO EXONS WITH PERIPHERAL INDELS)
            if indels >= allowed_indels and stretch_length < min_stretch_length:
                return homology

            # end of the exon
            if cds_seq[i] != gene_seq[i]:
                return homology

        # if loop is finished and homology not called -> default False        
        return homology


def check_synteny_5_prime(starting_position, locus_seq, gene_seq):
    """Checks synteny at 5'"""

    synteny = False
    UTR_length = 200
    stretch_length = 0
    matches = 0
    homology_threshold = 0.75

    # start from previous position
    starting_position -= 1

    # if contig too short, call lack of synteny
    if starting_position - UTR_length < 0:

        message = "\n..............Contig too short to check 5-prime synteny \n"
        print(message)
        logging.info(message)

        return synteny

    # if contig is long enough, estimate synteny
    for i in range(starting_position, starting_position - UTR_length + 1, -1):

        # skip indels
        if locus_seq[i] == "-" or gene_seq[i] == "-":
            continue

        # match
        elif locus_seq[i] == gene_seq[i]:
            matches += 1
            stretch_length += 1

        # mismatch
        elif locus_seq[i] != gene_seq[i]:
            stretch_length += 1

    # avoid divisions over 0
    if stretch_length == 0:

        message = "()()() stretch_length == 0"
        print(message)
        logging.info(message)

        return synteny

    # cannot call synteny if too many indels
    if stretch_length < 100:

        message = "\n..............Contig too short to check 5-prime synteny: only %s bp aligned \n" % str(stretch_length)
        print(message)
        logging.info(message)

        return synteny

    # check synteny
    if matches/stretch_length >= homology_threshold:
        synteny = True
        message = "\n..............5-prime syntenic: %s hamming distance\n" % str(matches/stretch_length)
        print(message)
        logging.info(message)

    else:
        synteny = False
        message = "\n..............5-prime NOT syntenic: %s hamming distance\n" % str(matches/stretch_length)
        print(message)
        logging.info(message)

    return synteny


def check_synteny_3_prime(starting_position, locus_seq, gene_seq):
    """Checks synteny at 3'"""

    synteny = False
    UTR_length = 200
    stretch_length = 0
    matches = 0
    homology_threshold = 0.75

    # if contig too short, call lack of synteny
    if starting_position + UTR_length > len(locus_seq):

        three_prime_synteny_message = "\n.............Alignment is too short to check 3-prime synteny \n"
        return synteny, three_prime_synteny_message

    # if contig is long enough, estimate synteny
    for i in range(starting_position, starting_position + UTR_length):

        # skip indels
        if locus_seq[i] == "-" or gene_seq[i] == "-":
            continue

        # match
        elif locus_seq[i] == gene_seq[i]:
            matches += 1
            stretch_length += 1

        # mismatch
        elif locus_seq[i] != gene_seq[i]:
            stretch_length += 1

    # avoid divisions over 0
    if stretch_length == 0:

        three_prime_synteny_message = "\n.............Contig too short to check 3-prime synteny: only %s bp aligned\n" \
                                                           % str(stretch_length)
        return synteny, three_prime_synteny_message

    # cannot call synteny if too many indels
    if stretch_length < 100:

        three_prime_synteny_message = "\n.............Contig too short to check 3-prime synteny: only %s bp aligned\n" \
                                      % str(stretch_length)
        return synteny, three_prime_synteny_message

    # check synteny
    if matches/stretch_length >= homology_threshold:
        synteny = True
        three_prime_synteny_message = "\n..............3-prime syntenic: %s hamming distance\n" \
                                      % str(matches/stretch_length)

    else:
        synteny = False
        three_prime_synteny_message = "\n..............3-prime NOT syntenic: %s hamming distance\n" \
                                      % str(matches/stretch_length)

    return synteny, three_prime_synteny_message


def check_retrotransposition(position, last_bp, cds_seq, locus_seq, gene_seq):  # runs only when intron is False
    """Checks for retrotransposition events"""

    retrotransposition = False

    if homology_check(position, last_bp, cds_seq, locus_seq, gene_seq) is True:
        retrotransposition = True
        return retrotransposition
    else:
        return retrotransposition


def reset_exon_parameters():
    """Resets parameters before the next exon is analysed"""

    exon = ""
    insertion = 0
    intron_at_5_prime = False
    intron_at_3_prime = False
    intron = False
    exon_missing = False
    last_bp = False
    big_insertion = False
    insertion_with_N = False

    return exon, insertion, intron_at_5_prime, intron_at_3_prime, \
         intron, exon_missing, last_bp, big_insertion, insertion_with_N


def check_insertion(contig_name, insertion, exon_number, ref_exons, insertion_with_N):
    """Checks for insertions"""

    big_insertion = False

    # do not allow Ns in insertions to avoid misalignment
    if insertion_with_N is True:
        big_insertion = True
        insertion_message = ">>>>> Contig %s exon %s contains %s bp insertion (Ns present -> not allowed) <<<<< \n" \
            % (contig_name, exon_number, str(insertion))
        print(insertion_message)
        logging.info(insertion_message)
        return big_insertion

    # allow insertions (even frameshifts) in last exons
    if exon_number == list(ref_exons.keys())[-1] and (insertion >= 180 or (insertion >= 10 and insertion % 3 != 0)):
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
