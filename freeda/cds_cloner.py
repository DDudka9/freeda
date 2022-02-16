#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 18:55:36 2021

@author: damian

Extracts single exons looking at exon/intron bounderies based on exons from ref species
Clones and stiches the exons into a final cds for given genome.

"""

from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from freeda import fasta_reader, pyinstaller_compatibility
import glob
import logging
import operator
import re
import shutil
import subprocess


def clone_cds(wdir, ref_species, preselected_exons_overhangs, most_intron_contigs, gene_name,
              genome_name, final_exon_number, ref_exons, MSA_path, aligner):
    """Clones cds for given genome based on the found exons"""
    
    # get a dictionary with all the contigs and how many exons they have
    all_contigs_dict = {}
    for contig in most_intron_contigs:
        contig_name, exon_nr = contig
        all_contigs_dict[contig_name] = exon_nr
    
    # estimate missing exons, duplicated exons and clone all possible cds
    missing_exons = []
    duplicated_exons = {}
    for e in range(1, len(preselected_exons_overhangs)+1):

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
        for contig in most_intron_contigs:
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
            # record all the analyzed contigs
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
                    contigs_hd_analyzed[winner] = "not analyzed"

            for winner in winners_hd_normalized:
                total_hd_normalized = 0
                exons = [exon_nr for exon_nr, cs in duplicated_exons.items() if winner in cs]
                for exon_nr in exons:
                    total_hd_normalized += hamming_distance_to_ref_species(wdir, ref_exons, exon_nr, winner,
                                                                           preselected_exons_overhangs, MSA_path,
                                                                           gene_name, aligner)
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
                
                message2 = "\n   ---->Contig " + winner + " was selected as most conserved"
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
        for contig in most_intron_contigs:
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
        
    # refill the list of most intron contigs
    most_intron_contigs_without_duplicates = []
    for contig in dup_exons_set:
        if contig in final_winners:
            c = contig, all_contigs_dict[contig]
            # append only the final winners that were never losers in any comparison
            most_intron_contigs_without_duplicates.append(c)
    # add to the list also all the contigs that didnt carry duplicated exons
    for contig in most_intron_contigs:
        c, exons_nr = contig
        if c not in dup_exons_set:
            most_intron_contigs_without_duplicates.append(contig)
            
    # finally clone cds using new most intron contig list
    most_intron_contigs = most_intron_contigs_without_duplicates
    cloned_cds = ""
    # make a scaffold for recording which exons was taken from which contig
    cds_composition = {exon: "" for exon, contigs in preselected_exons_overhangs.items()}
    
    # dictionary to collect sequences before finally stiching exons together
    pre_cloned_cds = {}
    # clone exon sequences from contigs (contigs with most intron exons are prioritized)
    for exon, contigs in preselected_exons_overhangs.items():
        # disregard missing exons
        if contigs != []:
            # look which contig has most intron exons
            for contig in most_intron_contigs:
                # look if a given contig contains a given exon
                for seq in preselected_exons_overhangs[exon]:
                    # if yes, and this exon wasnt yet cloned -> run MSA against ref species exon
                    if contig[0] == seq[0] and cds_composition[exon] == "":
                        in_filename, out_filename = generate_single_exon_MSA(wdir, ref_species, seq[1], contig[0], exon,
                                                                             gene_name, ref_exons, MSA_path, seq[2],
                                                                             aligner)
                        # collect the aligned sequences
                        aligned_seqs = collect_sequences(in_filename + out_filename)

                        # clone the locus_exon by:
                        # indexing positions in the alignement
                        ref_species_exon, locus_exon = index_positions_exons(aligned_seqs)

                        # converting the nucleotides into a string
                        locus_exon_string = "".join(locus_exon.values())
                        # mark that this exon was already cloned
                        cds_composition[exon] = contig[0]
                        pre_cloned_cds[exon] = contig[0], locus_exon_string
                        break
    
    # count how many exons are missing in the final cds cloned
    final_missing_exons_count = 0
    final_missing_exons = []
    
    # define which contigs were rejected
    final_exon_breakdown = {}

    # run again through cds_composition and draw final breakdown
    for exon, contig in cds_composition.items():

        if contig != "":
            final_exon_breakdown[exon] = contig
            continue

        else:
            final_exon_breakdown[exon] = "missing"
            final_missing_exons_count += 1
            final_missing_exons.append(exon)

    # clone the final cds
    cloned_cds = "".join([value[1] for key, value in pre_cloned_cds.items()])

    message1 = "\nCDS for %s gene in %s is composed of %s/%s exons from contigs: \n%s" \
            % (gene_name, genome_name, final_exon_number-final_missing_exons_count,
               final_exon_number, final_exon_breakdown)
    print(message1)
    logging.info(message1)
          
    message2 = "\nMissing exons: %s" % str(final_missing_exons).lstrip("[").rstrip("]")
    print(message2)
    logging.info(message2)
    
    return cloned_cds 


def hamming_distance_to_ref_species(wdir, ref_exons, exon_nr, winner, preselected_exons_overhangs,
                                    MSA_path, gene_name, aligner):
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
        f.write(">_" + gene_name + "gene\n")
        f.write(sequence[0][1])

    # prepare handles for msa
    in_filename = glob.glob(wdir + "/" + filename)[0]
    out_filename = filename.rstrip(".fasta") + "_aligned.fasta"

    # define which aligner is used
    if aligner == "mafft":

        cline = MafftCommandline(cmd=pyinstaller_compatibility.resource_path("mafft"),
                                 input=in_filename,
                                 thread=-1)  # thread -1 is suppose to automatically calculate physical cores
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

    # find hamming distance for them
    hd = hamming_distance_frameshift(ref_exon, compared_exon)  # compared_exon, exon_nr, ref_exons
    message = "   Hamming distance for contig " + winner + " exon nr: %s is %s" \
        % (str(exon_nr), str(hd))
    print(message)
    logging.info(message)
    
    # move both files to MSA_path
    shutil.move(in_filename, MSA_path + "Duplicated_exons/")
    shutil.move(out_filename, MSA_path + "Duplicated_exons/")

    #    return hd_normalized
    hd_normalized = hd / len(ref_exons[exon_nr][1])

    return hd_normalized
    
    
def index_positions_exons(seqs):
    """Indexes nucleotides for each exon into a dictionary"""

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
    """Checks if there is any frameshifts in the given exon. It cannot detect frameshifts at the end of contigs"""

    frameshift = False
    end_of_contig = False
    nr_of_insertions = 0
    nr_of_deletions = 0
    
    # check if there are dashes in ref_exon dictionary (insertion in compared_exon)
    if "-" in ref_exon.values():
        # count them
        count = 0
        lam = lambda pos: count + 1
        nr_of_insertions = len([lam(pos) for pos, value in ref_exon.items() if value == "-"])
        end_of_contig = True
        for nt in compared_exon.values():
            if nt == "-":
                end_of_contig = True
            else:
                end_of_contig = False
    
    if end_of_contig is True:
        frameshift = False
        return frameshift

    # check if there are dashes in compared_exon dictionary (deletion in compared_exon)
    if "-" in ref_exon.values() and end_of_contig is False:
        # count them
        count = 0
        lam = lambda pos: count + 1
        nr_of_deletions = len([lam(pos) for pos, value in compared_exon.items() if value == "-"])
    
    # make sure that indels altogether dont offset the frame
    if (nr_of_insertions - nr_of_deletions) % 3 != 0:
        frameshift = True
        return frameshift
    else:
        return frameshift


def hamming_distance_frameshift(p, q):  # exon_nr, ref_species_exons
    """Checks hamming distance between likely duplicated exons. Frameshifts are penalized,"""

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

    return len(mismatches)


def generate_single_exon_MSA(wdir, ref_species, seq, contig_name, exon_number, gene_name,
                             ref_exons, MSA_path, genomic_locus, aligner):
    """Writes into a file a reference single exon, part of the presumptive locus encoding it in a searched genome
    and a part of the reference genome encoding that exon"""
            
    filename = "exon_" + str(exon_number) + "_to_align.fasta"
    with open(filename, "w") as f:
        # write ref species exon first
        f.write(">_" + ref_species + "__" + gene_name + "_exon_" + str(exon_number))
        f.write("\n" + ref_exons[exon_number][1])
        # write the locus_exon second
        f.write("\n>_" + contig_name + "_" + gene_name + "_exon_" + str(exon_number))
        f.write("\n" + seq)
        # write the gene sequence third
        f.write("\n>_" + contig_name + "_gene")
        f.write("\n" + genomic_locus)
    in_filepath = MSA_path + "Single_exons_MSA/"
    # move the file for MSA
    shutil.move(wdir + filename, in_filepath)
    # run MSA and get the final filename of the aligned sequences
    out_filename = run_single_exon_msa(wdir, ref_species, in_filepath, str(exon_number), filename,
                                       aligner, ref_exons)

    return in_filepath, out_filename


def single_exon_mapping_checkpoint(wdir, ref_species, out_filename, in_filepath, exon_number, ref_exons):
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

    # eliminate misaligned exons by converting all bases to dashes (save original alignment as "_uncorrected")
    if 0.60 <= score <= 0.75:  # Cenpu Ay _LIPJ01001187.1__rev_Cenpu_exon_5 0.6363636363636364 -> legit

        # check if both flanks of the poorly aligned cloned exon are homologous to reference introns
        double_intron = check_intron_single_exon(ref_exon_aligned, cloned_exon_aligned, gene_aligned)

        if double_intron is False:

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


def check_intron_single_exon(ref_exon_aligned, cloned_exon_aligned, gene_aligned):
    """Checks intron of a single exon alignment. Exon overhangs are always 50bp long."""

    double_intron = True

    N_mismatches = 0
    N_indels = 0
    C_mismatches = 0
    C_indels = 0

    for position, bp in enumerate(cloned_exon_aligned):

        # check N-term intron
        if position < 50:

            # expecting intron regions (sometimes there are trailing true exonic bp that will also be penalized here)
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

        # check C-term intron
        if position > (len(cloned_exon_aligned) - 50):

            # expecting intron regions (sometimes there are trailing true exonic bp that will also be penalized here)
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

    if N_score < 0.75 or C_score < 0.75:
        double_intron = False

    return double_intron


def run_single_exon_msa(wdir, ref_species, in_filepath, exon_number, in_filename, aligner, ref_exons):
    """Runs MSA comparing reference species single exons and syntenic locus with presumptive exons"""

    out_filename = "aligned_" + "exon_" + str(exon_number) + ".fasta"

    # define which aligner is used
    if aligner == "mafft":

        cline = MafftCommandline(cmd=pyinstaller_compatibility.resource_path("mafft"),
                                 input=in_filepath + in_filename)
        # record standard output and standard error
        stdout, stderr = cline()
        # make a post-MSA file using out_filename
        with open(out_filename, "w") as f:
            f.write(stdout)
        # check for exon mapping errors -> eliminate such cloned exons
        single_exon_mapping_checkpoint(wdir, ref_species, out_filename, in_filepath, exon_number, ref_exons)
        # move the file to MSA_path
        shutil.move(out_filename, in_filepath)

    if aligner == "muscle":

        # outputs into working directory
        cmd = ['muscle', "-in", in_filepath + in_filename, "-quiet", "-maxiters", "2", "-out", out_filename]
        subprocess.call(cmd)
        # need to reorder seqs post msa
        fasta_reader.reorder_alignment(in_filepath + in_filename, wdir + out_filename)  # outputs to working dir

        # check for exon mapping errors -> eliminate such cloned exons
        single_exon_mapping_checkpoint(wdir, ref_species, out_filename, in_filepath, exon_number, ref_exons)
        shutil.move(out_filename, in_filepath)

    return out_filename


def collect_sequences(path):
    """Collects sequences from an alignment into a list."""

    # use biopython to read alignment file
    alignment = AlignIO.read(path, "fasta")
    # collect sequences in the alignment
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]

    return seqs

"""

    for position, bp in enumerate(ref_exon_aligned):

        if bp != "-" and bp == gene_aligned[position]:
            
            
        # allows for deletions and it is blind to insertions
        if bp != "-" \
                and cloned_exon_aligned[position] != ref_exon_aligned[position] \
                and cloned_exon_aligned[position] != "-":
            mapping_positions += 1

        # allows for deletions and it is blind to insertions
        elif bp != "-" and cloned_exon_aligned[position] == ref_exon_aligned[position]:
            matches += 1
            mapping_positions += 1


# penalizes indels
def single_exon_mapping_checkpoint(wdir, ref_species, out_filename, in_filepath, exon_number, ref_exons):
    Checks if single exons were properly mapped after exon finding

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
    exon_start = False
    exon_length = ref_exons[int(exon_number)][2]

    for position, bp in enumerate(ref_exon_aligned):

        if exon_length == 0:
            break

        # count bases in ref exon
        if bp != "-" and bp == gene_aligned[position]:
            exon_start = True
            exon_length -= 1

            # match
            if cloned_exon_aligned[position] == ref_exon_aligned[position]:
                matches += 1
                mapping_positions += 1

            # mismatch are penalized
            elif cloned_exon_aligned[position] != "-":
                mapping_positions += 1

            # deletions are moderately penalized
            elif cloned_exon_aligned[position] == "-":
                matches += 0.25
                mapping_positions += 1

        # insertion (exon start prevents counting insertions in intronic overhangs) - moderately penalized
        if exon_start and bp == "-" and cloned_exon_aligned[position] != "-":
            matches += 0.25
            mapping_positions += 1

        # misalignment -> reset matches
        if bp != "-" and bp != gene_aligned[position]:
            matches = 0

    score = matches/mapping_positions

    # eliminate misaligned exons by converting all bases to dashes (save original alignment as "uncorrected"
    if 0.60 <= score <= 0.75:  # Cenpu Ay _LIPJ01001187.1__rev_Cenpu_exon_5 0.6363636363636364 -> legit

        # check if both flanks of the poorly aligned cloned exon are homologous to reference introns
        double_intron = check_intron_single_exon(ref_exon_aligned, cloned_exon_aligned, gene_aligned)

        if double_intron is False:

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


"""