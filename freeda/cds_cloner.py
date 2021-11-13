#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 18:55:36 2021

@author: damian

Extracts single exons looking at exon/intron bounderies based on exons from ref species
Clones and stiches the exons into a final cds for given genome.

CURRENTLY (03_24_2021) THERE IS NO ALTERNATIVE STOP CODON SEARCH FUNCTION ->
# under assumption that Mm coding sequnce is would miss whatever insertion would
# happen in the other sequences and so PAML would be blind to it
# however this way a possible deletion in last exon would still be cloned
# because of "soft frameshift" rule that allows cloning frameshifts in last exon

"""

from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from freeda import fasta_reader
import glob
import logging
import operator
import re
import shutil
import subprocess


def clone_cds(wdir, ref_species, preselected_exons_overhangs, most_intronic_contigs, protein_name,
              genome_name, final_exon_number, ref_exons, MSA_path, aligner):
    
    # get a dictionary with all the contigs and how many exons they have
    all_contigs_dict = {}
    for contig in most_intronic_contigs:
        contig_name, exon_nr = contig
        all_contigs_dict[contig_name] = exon_nr
    
    # estimate missing exons, duplicated exons and clone all possible cds
    missing_exons = []
    duplicated_exons = {}
    for e in range(1, len(preselected_exons_overhangs)+1):
        # microexons are not expected in the preselected exons
        #if e in microexons:
        #    continue
        
        # if an exon is missing, add it to missing exons
        if preselected_exons_overhangs[e] == []:
            missing_exons.append(e)
            continue
        
        # if an exon is found in more than one contig
        if len(preselected_exons_overhangs[e]) > 1:
            contigs_with_duplicates = []
            # get contigs carrying duplicated exons
            for contig, duplicated_exon in enumerate(preselected_exons_overhangs[e]):
                contigs_with_duplicates.append(preselected_exons_overhangs[e][contig][0])
            duplicated_exons[e] = contigs_with_duplicates
    
    # get a set of all duplicated contigs
    dup_exons = [contigs for exon, contigs in duplicated_exons.items()]
    dup_exons_set = set([contig for l in dup_exons for contig in l])
    
    # find and preserve contigs with more exons if duplicated
    # this allows re-evaluation of duplicated contigs
    losing_contigs = set()
    winning_contigs = set()
    winning_and_losing_contigs = set()
    
    contigs_hd_analyzed = {}
    # for each duplicated exon, try to get a winner (most exons)
    for exon_nr, contigs in duplicated_exons.items():
        contigs_dict = {}
        for contig in most_intronic_contigs:
            contig_name, nr_of_exons = contig
            if contig_name in contigs:
                contigs_dict[contig_name] = nr_of_exons
        
        max_value = max(contigs_dict.values())
        winners = [contig_name for contig_name, nr_of_exons in contigs_dict.items() if nr_of_exons == max_value]
        # get the contig with most exons for that duplicated exon comparison
        winner = max(contigs_dict.items(), key=operator.itemgetter(1))[0]        
        
        # if there is only one winner
        if len(winners) == 1:
            winning_contigs.add(winner)
            # list the losers
            losers = [contig for contig in contigs if contig != winner]
            for loser in losers:
                losing_contigs.add(loser)
        
        # if any contig has been already analyzed for hd -> skip all
        elif any([contig_name for contig_name, nr_of_exons in contigs_dict.items() \
                if contig_name in contigs_hd_analyzed]):
            continue
        
        else:
            # record all the analysed contigs
            for contig_name, nr_exons in contigs_dict.items():
                contigs_hd_analyzed[contig_name] = 0
            message = "\nRunning hamming distance between contigs: " + str(winners)
            print(message)
            logging.info(message)
            
            # estimate total hd of each winner (normalized to exon length)
            winners_hd_normalized = {}
            for winner in winners:
                if contigs_dict[winner] == max_value:
                    winners_hd_normalized[winner] = 0
                else:
                    contigs_hd_analyzed[winner] = "not analysed"

            for winner in winners_hd_normalized:
                total_hd_normalized = 0
                exons = [exon_nr for exon_nr, cs in duplicated_exons.items() if winner in cs]
                for exon_nr in exons:
                    total_hd_normalized += hamming_distance_to_ref_species(wdir, ref_exons, exon_nr, winner,
                                                    preselected_exons_overhangs, MSA_path, protein_name, aligner)
                winners_hd_normalized[winner] = total_hd_normalized
                contigs_hd_analyzed[winner] = total_hd_normalized

            # get most conserved winner
            lowest_hd_normalized = min(winners_hd_normalized.values())
            winner = min(winners_hd_normalized.items(), key=operator.itemgetter(1))[0]
            for contig_name, nr_exons in contigs_dict.items():
                contigs_hd_analyzed[contig_name] = 0
            
            if lowest_hd_normalized == float("inf"):
                winner = ""
                message = "\n   ---->Contigs: " + str(winners) + "did NOT pass frameshift checkpoint"
                print(message)
                logging.info(message)
            
            else:
                message1 = "\nNormalized hamming distances are: " + str(winners_hd_normalized)
                print(message1)
                logging.info(message1)
                
                message2 = "\n   ---->Contig " + winner + \
                                " was selected as most conserved"
                print(message2)
                logging.info(message2)
            
            winning_contigs.add(winner)
            # list the losers
            losers = [contig for contig in contigs if contig != winner]
            for loser in losers:
                losing_contigs.add(loser)
                
    # list contigs that were both winners and losers
    for contig in dup_exons_set:
        if contig in winning_contigs and contig in losing_contigs:
            winning_and_losing_contigs.add(contig)
    
    # re-evaluate the duplicated exons 
    final_winners = set()
    for exon_nr, contigs in duplicated_exons.items():
        contigs_dict = {}
        for contig in most_intronic_contigs:
            contig_name, nr_of_exons = contig
            # skip contigs that are both winners and losers
            if contig_name in losing_contigs or contig_name in winning_and_losing_contigs:
                continue
            # add other contigs that should only be winners at this point
            if contig_name in contigs:
                contigs_dict[contig_name] = nr_of_exons
            # this contig doesnt have a duplicate -> skip
            else:
                continue
        # get the final winner for a given comparison
        try:
            winner = max(contigs_dict.items(), key=operator.itemgetter(1))[0]
            final_winners.add(winner)
        # if there is no contig in contig_dict (no winner for a given exon)
        except ValueError:
            continue
        
    # refill the list of most intronic contigs
    most_intronic_contigs_without_duplicates = []
    for contig in dup_exons_set:
        if contig in final_winners:
            c = contig, all_contigs_dict[contig]
            # append only the final winners that were never losers in any comparison
            most_intronic_contigs_without_duplicates.append(c)
    # add to the list also all the contigs that didnt carry duplicated exons
    for contig in most_intronic_contigs:
        c, exons_nr = contig
        if c not in dup_exons_set:
            most_intronic_contigs_without_duplicates.append(contig)
            
    # finally clone cds using new most intronic contig list
    most_intronic_contigs = most_intronic_contigs_without_duplicates
    cloned_cds = ""
    # make a scaffold for recording which exons was taken from which contig
    cds_composition = {exon: "" for exon, contigs in preselected_exons_overhangs.items()}
    
    # dictionary to collect sequences before finally stiching exons together
    pre_cloned_cds = {}
    exons_with_frameshifts = []
    # clone exon sequences from contigs (contigs with most intronic exons are prioritized)
    for exon, contigs in preselected_exons_overhangs.items():
        # disregard missing exons
        if contigs != []:
            # look which contig has most intronic exons
            for contig in most_intronic_contigs:
                # look if a given contig contains a given exon
                for seq in preselected_exons_overhangs[exon]:
                    # if yes, and this exon wasnt yet cloned -> run MSA against ref species exon
                    if contig[0] == seq[0] and cds_composition[exon] == "":
                        in_filename, out_filename = generate_single_exon_MSA(wdir, ref_species, seq[1], contig[0], exon,
                                                                    protein_name, ref_exons, MSA_path, seq[2], aligner)
                        # collect the aligned sequences
                        aligned_seqs = collect_sequences(in_filename + out_filename)

                        # clone the locus_exon by:
                        # indexing positions in the alignement
                        ref_species_exon, locus_exon = index_positions_exons(aligned_seqs)

                        # check frameshift (DISABLED)
                        frameshift = False
                        #frameshift = check_frameshift(ref_species_exon, locus_exon)
                        
                        # check if its not the last exon
                        if exon != list(ref_exons.keys())[-1] and frameshift is False:
                            # converting the nucleotides into a string
                            locus_exon_string = "".join(locus_exon.values())
                            # mark that this exon was already cloned
                            cds_composition[exon] = contig[0]
                            pre_cloned_cds[exon] = contig[0], locus_exon_string 
                            break 
                        
                        # if frameshift detected and not last exon
                        if frameshift is True and exon != list(ref_exons.keys())[-1]:
                            # log the frameshift
                            message = "   ...WARNING... : FRAMESHIFT detected in contig " \
                                    + contig[0] + " in exon nr: %s" % str(exon)
                            print(message)
                            logging.info(message)
                            # mark that this exon was already counted but not cloned
                            cds_composition[exon] = "FRAMESHIFT"
                            exons_with_frameshifts.append(exon)
                            pre_cloned_cds[exon] = contig[0], ""
                            break

                        # if last exon
                        if exon == list(ref_exons.keys())[-1]:
                            # check if there is a STOP codon
                            STOP = check_stop_codon(aligned_seqs[0][1], aligned_seqs[1][1], exon)
                            # converting the nucleotides into a string
                            locus_exon_string = "".join(locus_exon.values())
                            # mark that this exon was already cloned
                            cds_composition[exon] = contig[0]
                            pre_cloned_cds[exon] = contig[0], locus_exon_string
                            if frameshift is True:
                                exons_with_frameshifts.append(exon)
                                # log the frameshift
                                message = "   ...WARNING... : FRAMESHIFT detected in contig " \
                                    + contig[0] + " in the LAST exon nr: %s (allowed)" % str(exon)
                                print(message)
                                logging.info(message)
                                
                            break
    
    # count how many exons are missing in the final cds cloned
    final_missing_exons_count = 0
    final_missing_exons = []
    
    # define which contigs were rejected
    rejected_contigs = set()
    final_exon_breakdown = {}
    for exon, contig in cds_composition.items():
        if contig == "FRAMESHIFT":
            rejected_contigs.add(pre_cloned_cds[exon][0])
    
    # run again through cds_composition and draw final breakdown
    for exon, contig in cds_composition.items():
        if contig in rejected_contigs:
            final_exon_breakdown[exon] = pre_cloned_cds[exon][0] + \
                " contains FRAMESHIFT in exons: " + str(exons_with_frameshifts) \
                                                    + " (not cloned)"
            continue
        if "FRAMESHIFT" not in contig and contig != "":
            final_exon_breakdown[exon] = contig
            continue
        else:
            final_exon_breakdown[exon] = "missing"
            final_missing_exons_count += 1
            final_missing_exons.append(exon)

    # clone the final cds
    cloned_cds = "".join([value[1] for key, value in pre_cloned_cds.items()])

    message1 = "\nCDS for %s protein in %s is composed of %s/%s exons from contigs: \n%s" \
            % (protein_name, genome_name, final_exon_number-final_missing_exons_count,
               final_exon_number, final_exon_breakdown)
    print(message1)
    logging.info(message1)
          
    message2 = "\nMissing exons: %s" % str(final_missing_exons).lstrip("[").rstrip("]")
    print(message2)
    logging.info(message2)
    
    #message3 = "\nDuplicated exons (contigs with RETRO_score >= 0.4 were excluded): %s" \
    #    % str(duplicated_exons.keys()).lstrip("dict_keys([").rstrip("])")
    #print(message3)
    #logging.info(message3)
    
    return cloned_cds 


def hamming_distance_to_ref_species(wdir, ref_exons, exon_nr, winner, preselected_exons_overhangs,
                                    MSA_path, protein_name, aligner):
    """Compares ref exon to the presumptive exon in order to pick between >1 exon candidate -> picks more conserved."""

    # determine the name for the file to MSA
    filename = "duplicated_exon_" + str(exon_nr) + "_in_" + winner + "_hamming_distance.fasta"

    # prepare a file to run MSA
    with open(filename, "w") as f:
        # write the reference species exon first
        f.write(ref_exons[exon_nr][0] + "\n")
        f.write(ref_exons[exon_nr][1] + "\n")
        # for each winner write in an exon
        e = preselected_exons_overhangs[exon_nr]
        sequence = [(seq, genomic) for (contig_name, seq, genomic) in e if contig_name == winner]
        f.write(">_" + "exon_" + str(exon_nr) + "_" + str(winner) + "\n")
        f.write(sequence[0][0] + "\n")
        f.write(">_" + protein_name + "gene\n")
        f.write(sequence[0][1])

    # prepare handles for msa
    in_filename = glob.glob(wdir + "/" + filename)[0]
    out_filename = filename.rstrip(".fasta") + "_aligned.fasta"

    # define which aligner is used
    if aligner == "mafft":

        cline = MafftCommandline(input=in_filename, thread=-1)  # thread -1 is suppose to automatically
                                                                # calculate physical cores
        # record standard output and standard error
        stdout, stderr = cline()
        # make a post-MSA file using out_filename
        with open(out_filename, "w") as f:
            f.write(stdout)

    if aligner == "muscle":

        cmd = ['muscle', "-in", in_filename, "-quiet", "-maxiters", "2", "-out", out_filename]
        subprocess.call(cmd)
        # need to reorder seqs post msa
        fasta_reader.reorder_alignment(in_filename, out_filename)

    path = wdir + "/" + out_filename
    # use biopython aligment reader module
    alignment = AlignIO.read(path, "fasta")
    # collect sequences in the alignment
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]
    # get trimmed and indexed versions of the reference species exon and compared exon
    ref_exon, compared_exon = index_positions_exons(seqs)
    # check if there is no frame shift in the given exon
    frameshift = check_frameshift(ref_exon, compared_exon)
    if frameshift is True:
        message = "   ...WARNING... : FRAMESHIFT detected in contig " + winner + \
            " in exon nr: %s" % str(exon_nr)
        print(message)
        logging.info(message)
    
    
    # write a checkpoint that will limit exon checking to only overlapping exons in winner contigs
    
    
    # find hamming distance for them
    hd = hamming_distance_frameshift(ref_exon, compared_exon, exon_nr, ref_exons)
    message = "   Hamming distance for contig " + winner + " exon nr: %s is %s" \
        % (str(exon_nr), str(hd))
    print(message)
    logging.info(message)
    
    # move both files to MSA_path
    shutil.move(in_filename, MSA_path + "Duplicated_exons/")
    shutil.move(out_filename, MSA_path + "Duplicated_exons/")

    # if there is a frameshift and its not the last exon, reject this contig
    if frameshift is True and exon_nr != list(ref_exons)[-1]:
        hd_normalized = float("inf")
        
        return hd_normalized
    else:
        # normalize hd to ref species exon (prior MSA; without indels)
        hd_normalized = hd/len(ref_exons[exon_nr][1])
        return hd_normalized    
    
    
def index_positions_exons(seqs):

    ref_exon = seqs[0][1]
    compared_exon = seqs[1][1]
    # get position of the last dash before the start of the ref species exon
    first_nucl = re.search(r"[ACTG]", seqs[0][1]).span()[0]
    # get position of the last character in ref species exon
    last_nucl = max([len(seqs[0][1].rsplit(n, 1)[0]) for n in "ACTG"])
    # trim ref exon (includes indels)
    trimmed_ref_exon = ref_exon[first_nucl:last_nucl+1]
    trimmed_compared_exon = compared_exon[first_nucl:last_nucl+1]
    # get indexes for trimmed exons
    trimmed_ref_exon_indexed = {}
    trimmed_compared_exon_indexed = {}
    for n in range(len(trimmed_ref_exon)):
        trimmed_ref_exon_indexed[n] = trimmed_ref_exon[n]
        trimmed_compared_exon_indexed[n] = trimmed_compared_exon[n]
    
    return trimmed_ref_exon_indexed, trimmed_compared_exon_indexed


def check_frameshift(ref_exon, compared_exon):

    frameshift = False
    end_of_contig = False
    oe = ref_exon
    ce = compared_exon
    nr_of_insertions = 0
    nr_of_deletions = 0
    
    # check if there are dashes in oe dictionary (insertion in ce)
    if "-" in oe.values():
        # count them
        count = 0
        lam = lambda pos: count + 1
        nr_of_insertions = len([lam(pos) for pos, value in oe.items() if value == "-"])
        
        # check if the dashes are in the end of the exon (likely end of contig)
        
        # THIS FUNCTION CANNOT DETECT FRAMESHIFTS AT THE END OF CONTIGS (ex. Cdk5rap2 An exon 11 NC_047662.1)
        # track back last 3 positions -> trim if frameshift detected -> clone the rest
        
        end_of_contig = True
        for nt in ce.values():
            if nt == "-":
                end_of_contig = True
            else:
                end_of_contig = False
    
    if end_of_contig == True:
        frameshift = False
        return frameshift

    # check if there are dashes in ce dictionary (deletion in ce)
    if "-" in oe.values() and end_of_contig == False:
        # count them
        count = 0
        lam = lambda pos: count + 1
        nr_of_deletions = len([lam(pos) for pos, value in ce.items() if value == "-"])
    
    # make sure that indels altogether dont offset the frame
    if (nr_of_insertions - nr_of_deletions) % 3 != 0:
        frameshift = True
        return frameshift
    else:
        return frameshift


def hamming_distance_frameshift(p, q, exon_nr, ref_species_exons):

    indels_p = 0
    indels_q = 0
    mismatches = []
    for i in range (len(p)):
        if p[i] != q[i]:
            mismatches.insert(0,i)
        if p[i] == "-":
            indels_p = indels_p + 1
        if q[i] == "-":
            indels_q = indels_q + 1
    
    # check if the dashes are in the end of the exon (likely end of contig)
    end_of_contig = True
    for nt in q.values():
        if nt == "-":
            end_of_contig = True
        else:
            end_of_contig = False
    
    # trigger highest possible hd if frameshift detected        
    if end_of_contig == False:
        if (indels_p % 3 != 0 or indels_q % 3 != 0) and exon_nr != list(ref_species_exons.keys())[-1]:
            
            return float("inf")

    return len(mismatches)


def generate_single_exon_MSA(wdir, ref_species, seq, contig_name, exon_number, protein_name,
                             ref_exons, MSA_path, genomic_locus, aligner):
    """Writes into a file a reference single exon, part of the presumptive locus encoding it in a searched genome
    and a part of the reference genome encoding that exon"""
            
    filename = "exon_" + str(exon_number) + "_to_align.fasta"
    with open(filename, "w") as f:
        # write ref species exon first
        f.write(">_" + ref_species + "__" + protein_name + "_exon_" + str(exon_number))
        f.write("\n" + ref_exons[exon_number][1])
        # write the locus_exon second
        f.write("\n>_" + contig_name + "_" + protein_name + "_exon_" + str(exon_number))
        f.write("\n" + seq)
        # write the gene sequence third
        f.write("\n>_" + contig_name + "_gene")
        f.write("\n" + genomic_locus)
    in_filepath = MSA_path + "Single_exons_MSA/"
    # move the file for MSA
    shutil.move(wdir + filename, in_filepath)
    # run MSA and get the final filename of the aligned sequences
    out_filename = run_single_exon_msa(wdir, ref_species, in_filepath, str(exon_number), filename, aligner)

    return in_filepath, out_filename


def single_exon_mapping_checkpoint(wdir, ref_species, out_filename, in_filepath):
    """Checks if single exons were properly mapped after exon finding."""

    # get dictionary from single exon alignment
    all_seq_dict = fasta_reader.alignment_file_to_dict(wdir, ref_species, out_filename)

    # write ref exon and cloned exon into str variables (with dashes)
    for i, header in enumerate(all_seq_dict):
        if i == 0:
            ref_exon_aligned = all_seq_dict[header]
        elif i == 1:
            cloned_exon_header = header
            cloned_exon_aligned = all_seq_dict[cloned_exon_header]
        else:
            gene_header = header
            gene_aligned = all_seq_dict[gene_header]

    # check for unreasonable insertions
    if len(cloned_exon_aligned.replace("-", "")) / len(gene_aligned.replace("-", "")) > 2:
        message = "\n...WARNING... : Exon : %s : Contains a huge insertions " \
                  "-> eliminated" % cloned_exon_header
        print(message)
        logging.info(message)
        uncorrected_filename = "_uncorrected_" + out_filename
        shutil.copy(wdir + out_filename, wdir + uncorrected_filename)
        shutil.move(wdir + uncorrected_filename, in_filepath)
        cloned_exon = "-" * len(cloned_exon_aligned)
        all_seq_dict[cloned_exon_header] = cloned_exon

        # write dict back to file
        with open(out_filename, "w") as f:
            for header, seq in all_seq_dict.items():
                f.write(header + "\n")
                f.write(seq + "\n")
            return

    # score mismatches
    matches = 0
    mapping_positions = 0
    for position, bp in enumerate(ref_exon_aligned):

        # allows for deletions and it is blind to insertions
        if bp != "-" \
                and cloned_exon_aligned[position] != ref_exon_aligned[position] \
                and cloned_exon_aligned[position] != "-":
            mapping_positions += 1

        elif bp != "-" and cloned_exon_aligned[position] == ref_exon_aligned[position]:
            matches += 1
            mapping_positions += 1

    score = matches/mapping_positions

    # eliminate misaligned exons by converting all bases to dashes (save original alignment as "uncorrected"
    if 0.60 <= score <= 0.75:  # Cenpu Ay _LIPJ01001187.1__rev_Cenpu_exon_5 0.6363636363636364 -> legit

        # check if both flanks of the poorly aligned cloned exon are homologous to reference introns
        double_introny = check_introny_single_exon(ref_exon_aligned, cloned_exon_aligned, gene_aligned)

        if double_introny is False:

            message = "\n...WARNING... : Exon : %s : Questionable alignment score : %s (good > 0.75, questionable " \
                      "= 0.60-0.75, poor < 0.60) and lack of two flanking introns " \
                      "-> eliminated" % (cloned_exon_header, score)
            print(message)
            logging.info(message)
            uncorrected_filename = "_uncorrected_" + out_filename
            shutil.copy(wdir + out_filename, wdir + uncorrected_filename)
            shutil.move(wdir + uncorrected_filename, in_filepath)
            cloned_exon = "-" * len(cloned_exon_aligned)
            all_seq_dict[cloned_exon_header] = cloned_exon

            # write dict back to file
            with open(out_filename, "w") as f:
                for header, seq in all_seq_dict.items():
                    f.write(header + "\n")
                    f.write(seq + "\n")
            return

        # low score but introns detected at both ends
        else:
            message = "\n...WARNING... : Exon : %s : Questionable alignment score : %s (good > 0.75, questionable " \
                      "= 0.60-0.75, poor < 0.60) but flanked by two syntenic introns " \
                      "-> acceptable" % (cloned_exon_header, score)
            print(message)
            logging.info(message)
            return

    if score < 0.60:

        message = "\n...WARNING... : Exon : %s : Poor alignment score : %s " \
                  "(good > 0.75, questionable = 0.60-0.75, poor < 0.60)" \
                  " -> eliminated" % (cloned_exon_header, score)
        print(message)
        logging.info(message)
        uncorrected_filename = "_uncorrected_" + out_filename
        shutil.copy(wdir + out_filename, wdir + uncorrected_filename)
        shutil.move(wdir + uncorrected_filename, in_filepath)
        cloned_exon = "-" * len(cloned_exon_aligned)
        all_seq_dict[cloned_exon_header] = cloned_exon

        # write dict back to file
        with open(out_filename, "w") as f:
            for header, seq in all_seq_dict.items():
                f.write(header + "\n")
                f.write(seq + "\n")
        return

    if score > 0.75:
        message = "\nExon : %s Good alignment score : %s (good > 0.75, questionable = 0.60-0.75, poor < 0.60) " \
                  "-> accepted" % (cloned_exon_header, score)
        print(message)
        logging.info(message)
        return


def check_introny_single_exon(ref_exon_aligned, cloned_exon_aligned, gene_aligned):
    """Checks introny of a single exon alignment. Exon overhangs are always 50bp long."""

    double_introny = True

    N_mismatches = 0
    N_indels = 0
    C_mismatches = 0
    C_indels = 0

    for position, bp in enumerate(cloned_exon_aligned):

        # check N-term introny
        if position < 50:

            # expecting intronic regions (sometimes there are trailing true exonic bp that will also be penalized here)
            if ref_exon_aligned[position] != "-":
                N_mismatches += 1

            # allow max 10 indels
            elif N_indels < 10:

                # deletion
                if bp == "-" and gene_aligned[position] != "-":
                    N_indels += 1
                # insertion
                elif bp != "-" and gene_aligned[position] == "-":
                    N_indels += 1
                # mismatch
                elif bp != gene_aligned[position]:
                    N_mismatches += 1
                # match
                else:
                    pass

            # past indel margin
            elif N_indels >= 10:

                # deletion
                if bp == "-":
                    N_mismatches += 1
                # insertion
                elif bp != gene_aligned[position]:
                    N_mismatches += 1
                # match
                else:
                    pass

        # check C-term introny
        if position > (len(cloned_exon_aligned) - 50):

            # expecting intronic regions (sometimes there are trailing true exonic bp that will also be penalized here)
            if ref_exon_aligned[position] != "-":
                C_mismatches += 1

            # allow max 10 indels
            elif C_indels < 10:

                # deletion
                if bp == "-" and gene_aligned[position] != "-":
                    C_indels += 1
                # insertion
                elif bp != "-" and gene_aligned[position] == "-":
                    C_indels += 1
                # mismatch
                elif bp != gene_aligned[position]:
                    C_mismatches += 1
                # match
                else:
                    pass

            # past indel margin
            elif C_indels > 10:

                # deletion
                if bp == "-":
                    C_mismatches += 1
                # insertion
                elif bp != gene_aligned[position]:
                    C_mismatches += 1
                # match
                else:
                    pass

    N_score = (50 - N_mismatches) / 50
    message = "\n    N_score : %s " % N_score
    print(message)
    logging.info(message)

    C_score = (50 - C_mismatches) / 50
    message = "    C_score : %s " % C_score
    print(message)
    logging.info(message)

    if (N_score or C_score) < 0.75:
        double_introny = False

    return double_introny


def run_single_exon_msa(wdir, ref_species, in_filepath, exon_number, in_filename, aligner):
    """Runs MSA comparing reference species single exons and syntenic locus with presumptive exons"""

    out_filename = "aligned_" + "exon_" + str(exon_number) + ".fasta"

    # define which aligner is used
    if aligner == "mafft":

        cline = MafftCommandline(input=in_filepath + in_filename)
        # record standard output and standard error
        stdout, stderr = cline()
        # make a post-MSA file using out_filename
        with open(out_filename, "w") as f:
            f.write(stdout)
        # check for exon mapping errors -> eliminate such cloned exons
        single_exon_mapping_checkpoint(wdir, ref_species, out_filename, in_filepath)
        # move the file to MSA_path
        shutil.move(out_filename, in_filepath)

    if aligner == "muscle":

        # outputs into working directory
        cmd = ['muscle', "-in", in_filepath + in_filename, "-quiet", "-maxiters", "2", "-out", out_filename]
        subprocess.call(cmd)
        # need to reorder seqs post msa
        fasta_reader.reorder_alignment(in_filepath + in_filename, wdir + out_filename)  # outputs to working dir

        # check for exon mapping errors -> eliminate such cloned exons
        single_exon_mapping_checkpoint(wdir, ref_species, out_filename, in_filepath)
        shutil.move(out_filename, in_filepath)

    return out_filename


def collect_sequences(path):

    # use biopython to read alignment file
    alignment = AlignIO.read(path, "fasta")
    # collect sequences in the alignment
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]

    return seqs


def check_stop_codon(last_ref_species_exon, last_locus_exon, exon_number): # works well

    # define the STOP codon based on the last ref species exon
    TAG_codon_pos = last_ref_species_exon.rfind("TAG")
    TGA_codon_pos = last_ref_species_exon.rfind("TGA")
    TAA_codon_pos = last_ref_species_exon.rfind("TAA")
    STOP_codon_pos = max(TAG_codon_pos, TGA_codon_pos, TAA_codon_pos)
    STOP = False
    STOP_locus_exon = last_locus_exon[STOP_codon_pos : STOP_codon_pos + 3]

    # check if the STOP codon is present
    if STOP_locus_exon == "TAG" or STOP_locus_exon == "TGA" or STOP_locus_exon == "TAA":
        STOP = True
        message = ("\nSTOP codon detected in LAST exon (%s)" % exon_number)
        print(message)
        logging.info(message)
    else:
        message = ("\nSTOP codon NOT detected in LAST exon (%s)" % exon_number)
        print(message)
        logging.info(message)

    return STOP


"""
        # check C-term introny
        if position > (len(cloned_exon_aligned) - 50):
            # match
            if bp == gene_aligned[position]:
                C_matches += 1
            # deletion in locus
            elif bp == "-" and C_indels < 10:
                C_indels += 1
            # insertion in locus
            elif bp != "-" and C_indels < 10:
                C_indels += 1

    N_score = N_matches / (50 - N_indels)
    print("    N_score : %s " % N_score)
    C_score = C_matches / (50 - C_indels)
    print("    C_score : %s " % C_score)
    
"""