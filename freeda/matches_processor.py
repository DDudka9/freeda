#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:40:58 2021

@author: damian

Extracts fasta sequence based on the matches dataframe and writes it into 
a fasta file together with cds and genomic locus of the protein from ref species.
Ready to MSA.

"""

from Bio.SeqRecord import SeqRecord
import logging
import shutil
import os
import pybedtools
import glob


def process_matches(ref_species, wdir, matches, cds_seq, gene_seq, result_path, gene, genome_name, genome_index,
                    all_genes_dict=None):
    """Processes blast matches and prepares them for the alignment with reference cds and gene"""

    # make a no-duplicates list of contig names
    sseqids = set(matches["sseqid"].tolist())
    # for a give contig:
    for contig in sseqids:
        # get exons as dataframe
        base_name = contig.split("__")[0]
        exons, bed_name, fasta_name = get_exons(base_name, matches, contig)
        # check for exons detected on reversed strand to pseudohap
        validated_exons, nr_of_exons_reversed, start, end, rev = check_strand(exons)
        # log how many reversed exons were detected
        log_strand(validated_exons, nr_of_exons_reversed, contig)
        # sort exons
        sorted_exons = sort_exons(validated_exons, rev)
        # expand exons
        expanded_exons = expand_exons(sorted_exons)
        # make a bed file with expanded exons or not expanded exons if cannot expand
        bed_object, expanded_bed_object = make_bed_file(sorted_exons, expanded_exons,
                                                        bed_name, result_path, gene, genome_name)
        # convert the bed file into fasta file
        make_fasta_file(wdir, fasta_name, bed_object, expanded_bed_object, genome_name)
        fasta_path = process_fasta_file(fasta_name, contig, result_path, gene, genome_name, rev)
        get_contig_locus(ref_species, contig, gene, genome_name, fasta_path, start, end, genome_index, rev,
                         all_genes_dict)
        MSA_path = generate_files_to_MSA(contig, cds_seq, gene_seq, fasta_path)

    message = "\nAnalysing all contigs ..."
    print(message)
    logging.info(message)

    return MSA_path


def get_exons(base_name, matches, contig):
    """Groups matches in presumable exons"""

    # group all entries from the dataframe for this contig
    grouped = matches.groupby(matches.sseqid)
    # make contig name a string
    contig = str(contig)
    message = "Matches found on contig : " + contig
    print(message)
    logging.info(message)
    # make a dataframe from these entries
    exons = grouped.get_group(contig)
    # bring back ref name of the contig for each match in exons (base_name)
    for i in exons.index:
        exons.at[i, "sseqid"] = base_name
    # prepare bed and fasta file names and path directories for the contig:
    bed_name = "_" + contig + ".bed"
    fasta_name = "_" + contig + ".fasta"
    # write dataframe to a file
    exons.to_csv(bed_name, sep=("\t"), index=False)
    # check if match is on reverse strand; swap "sstart" and "send" if true:
    # THIS DOESNT SEEM TO DO MUCH - DO I NEED THESE TWO LINES? :
    default = ["sseqid", "sstart", "send", "qseqid"]
    exons = exons.reindex(columns=default)

    return exons, bed_name, fasta_name


def check_strand(exons):
    """Checks strands of the presumable exons and reindexes the dataframe accordingly"""

    rev = False
    nr_of_exons_reversed = 0
    # select only "sstart" and "send" columns from exons dataframe
    columns = ["sstart", "send"]
    sstart_send = exons[columns]
    # find min value among all integers in "sstart_send" dataframe
    start = sstart_send.to_numpy().min()
    # find max value among all integers in "sstart_send" dataframe
    end = sstart_send.to_numpy().max()
    # iterate over rows in dataframe
    for index, row in exons.iterrows():
        # check for opposite strand match
        if row["sstart"] > row["send"]:
            # global rev
            rev = True
            nr_of_exons_reversed += 1
            # reindex columns in this row
            reverse = ["sseqid", "send", "sstart", "qseqid"]
            row2 = row.reindex(reverse)
            # rename columns in this row
            row2 = row2.rename({"send": "sstart", "sstart": "send"})
            # add that row to exons dataframe
            exons = exons.append(row2)
            # ckeck for duplicated indexes
            if exons.index.has_duplicates:
                # keep the first duplicated index (the reversed one)
                exons = exons[exons.index.duplicated(keep="first")]

    return exons, nr_of_exons_reversed, start, end, rev


def log_strand(exons, rows_reversed, contig):
    """Checkpoint for contigs with matches on opposite strands - which would be incorrect"""

    # check if there are matches on both strands in a contig
    if len(exons) > rows_reversed and rows_reversed != 0:
        message = "contig: " + contig + " has matches on both strands!"
        # log and print warning
        print(message)
        logging.info(message)


def sort_exons(exons, rev):
    """Sorts presumable exons by start position"""

    # check for opposite strand match
    if rev:
        # sort matches descending
        sorted_exons = exons.sort_values(by=["sstart"], ascending=True)
    else:
        # sort matches ascending
        sorted_exons = exons.sort_values(by=["sstart"], ascending=False)
    return sorted_exons


def expand_exons(exons):
    """Expands presumable exons as blast often misses some bp"""

    # need to make an idependent copy of dataframe "exons"
    expanded_exons = exons.copy()
    # generate Series for "sstart" subtracting 3bp (avoid blast cutting a codon)
    extended_sstart = expanded_exons["sstart"].sub(3)
    # add this Series to Dataframe
    expanded_exons.insert(1, column="extended_sstart", value=extended_sstart)
    # remove old "sstart" Series from Dataframe
    expanded_exons = expanded_exons.drop(columns="sstart")
    # rename "extended_sstart" Series to "sstart"
    expanded_exons = expanded_exons.rename(columns={"extended_sstart": "sstart"})
    # generate Series for "send" adding 3bp (avoid blast cutting a codon)
    extended_send = expanded_exons["send"].add(3)
    # add this Series to Dataframe
    expanded_exons.insert(2, column="extended_send", value=extended_send)
    # remove old "send" Series from Dataframe
    expanded_exons = expanded_exons.drop(columns="send")
    # rename "extended_send" Series to "send"
    expanded_exons = expanded_exons.rename(columns={"extended_send": "send"})

    return expanded_exons


def make_bed_file(exons, expanded_exons, bed_name, result_path, gene, genome_name):
    """Makes a bed file with coordinates for each contig (containing presumable exons)"""

    # make bed file
    expanded_bed_object = pybedtools.BedTool(bed_name)
    expanded_bed_object = pybedtools.bedtool.BedTool.from_dataframe(expanded_exons, header=False, index=False)
    # in case match is at end/beginning of contig dont add bp
    bed_object = pybedtools.bedtool.BedTool.from_dataframe(exons, header=False, index=False)
    # move bed file to Bed folder
    bed_path = result_path + gene + "/" + genome_name + "/Bed/above_threshold"
    shutil.move(bed_name, bed_path)

    return bed_object, expanded_bed_object


def make_fasta_file(wdir, fasta_name, bed_object, expanded_bed_object, genome_name):
    """Makes a fasta file for each contig (containing presumable exons) from a given genome based on bed file"""

    genome_dir = wdir + "Genomes/"

    # get fasta file from genome
    try:
        fasta_file = expanded_bed_object.sequence(fi=genome_dir + genome_name + ".fasta")
        fasta_file = fasta_file.save_seqs(fasta_name)

    # assuming that exception comes from match at end/beginning of contig
    except:
        fasta_file = bed_object.sequence(fi=genome_dir + genome_name + ".fasta")
        fasta_file = fasta_file.save_seqs(fasta_name)


def process_fasta_file(fasta_name, contig, result_path, gene, genome_name, rev):
    """Stiches presumable exons together (reverse complements if needed) and writes back into a fasta file"""

    # reverse complement if "sstart" > "send" (marked as rev == True)
    fasta_path = result_path + gene + "/" + genome_name + "/Fasta/above_threshold"

    if rev:
        # reverse complements and stiches presumable exons together
        rev_comp_filename = reverse_complement(fasta_name, gene, genome_name, contig, rev)
        # move the reverse complemented fasta file to Fasta folder
        shutil.move(rev_comp_filename, fasta_path)
        os.remove(fasta_name)

    else:
        # stich exons without reverse complementing
        exons_stich(fasta_name, gene, genome_name, rev=False)
        # move fasta file to Fasta folder
        shutil.move(fasta_name, fasta_path)

    return fasta_path


def exons_stich(fasta_name, gene, genome_name, rev=False):
    """Stiches presumable exons together and writes into a file"""

    header, seq = read_seq(fasta_name, rev)
    write_seq(header, seq, gene, genome_name, fasta_name)


def read_seq(fasta_name, rev=False):
    """Reads a fasta file with presumbale exon sequence"""

    if not rev:

        with open(fasta_name) as i:
            seq = ""
            head = ""
            for line in i.readlines():
                if line.startswith(">"):
                    head = line.lstrip(">").rstrip("\n") + head
                else:
                    seq = line.rstrip("\n").upper() + seq

        header = ">_" + head

        return header, seq

    else:

        with open(fasta_name) as i:
            seq = ""
            header = ">"
            for line in i.readlines():
                if line.startswith(">"):
                    head = line.lstrip(">").rstrip("\n")
                    header = header + "_" + head
                else:
                    seq = seq + line.rstrip("\n").upper()

        return header, seq


def write_seq_rev_comp(header, complement, gene, genome_name, contig):
    """Writes sequence of the reverse compemented contig into a fasta file"""

    # check and log if "N" bases are present
    non_ACGT = ["N", "Y", "R", "W", "S", "K", "M", "D", "H", "V", "B", "X"]

    if [n for n in complement if n in non_ACGT] != []:
        side_note = "      ...WARNING... : non-ACGT bases in the reverse complemented sequence."
        print(side_note)
        logging.info(side_note)

    o = open("Rev_comp_" + contig + ".fasta", "w")
    # generate a global filename for this sequence; start from default empty string
    # by converting file name to string that "shutil.move()" function can use
    rev_comp_filename = str(o.name)
    h = ">_Rev_comp_" + gene + "_" + genome_name + header.lstrip(">")
    o.write(h.rstrip("\n"))
    o.write("\n" + complement)
    o.close()

    return rev_comp_filename


def write_seq(header, seq, gene, genome_name, fasta_name):
    """Writes sequence of the contig into a fasta file"""

    # check and log if "N" bases are present
    non_ACGT = ["N", "Y", "R", "W", "S", "K", "M", "D", "H", "V", "B", "X"]

    if [n for n in seq if n in non_ACGT] != []:
        side_note = "    ...WARNING... : non-ACGT bases in the sequence."
        print(side_note)
        logging.info(side_note)

    o = open(fasta_name, "w")
    h = ">_" + gene + "_" + genome_name + header.lstrip(">")
    o.write(h.rstrip("\n"))
    o.write("\n" + seq)
    o.close()


def reverse_complement(fasta_name, gene, genome_name, contig, rev):
    """Makes a reverse complement"""

    rules = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
             "Y": "R", "R": "Y", "W": "W", "S": "S", "K": "M", "M": "K",
             "D": "H", "H": "D", "V": "B", "B": "V", "X": "X"}

    header, seq = read_seq(fasta_name, rev)
    complement = ""
    for base in seq:
        complement = rules[base] + complement

    rev_comp_filename = write_seq_rev_comp(header, complement, gene, genome_name, contig)

    return rev_comp_filename


def get_contig_locus(ref_species, contig, gene, genome_name, fasta_path, start, end, genome_index, rev,
                     all_genes_dict=None):
    """Generates a fasta file for given contig locus"""

    # naming of the file and header will depend on if contig was rev_comp or not
    if not rev:
        file_name = "Contig_locus_" + contig + ".fasta"
        header = ">_" + file_name.rstrip(".fasta") + "_" + gene + "_" + genome_name
    else:
        file_name = "Rev_comp_contig_locus_" + contig + ".fasta"
        header = ">_" + file_name.rstrip(".fasta") + "_" + gene + "_" + genome_name
    # get Seq object out of the indexed genome using contig string as a key
    # pull the contig by its base name
    base_name = contig.split("__")[0]
    seq = genome_index[str(base_name)].seq
    prefix, suffix = get_prefix_suffix(ref_species, start, len(seq), gene, all_genes_dict)
    # try trim contig and add 10000bp overhangs for later introny check
    seq = seq[start - prefix: end + suffix]
    # if full contig lacks enough bp on either end, the Seq object
    # will have length 0 and needs to be created again from genome.index
    # if seq1 len=0 then seq2 surely too, no need to check

    if rev:
        # if contig is too short, proceed without trimming and overhangs
        # reverse_complement the Seq object
        seq = seq.reverse_complement()
    # transform Seq object into SeqRecord object for the sake of format method
    record = SeqRecord(seq)
    # use format method of SeqRecord to get sequence string in fasta format
    sequence = record.format("fasta")
    # generate a file for full contig
    o = open(file_name, "w")
    o.write(header.rstrip("><unknown id> <unknown description>\n"))
    o.write("\n" + sequence.lstrip("><unknown id> <unknown description>\n"))
    o.close()
    # move file to fasta_path
    shutil.move(file_name, fasta_path)


def get_prefix_suffix(ref_species, start, seq_length, gene, all_genes_dict):
    """Generates extension for each contig"""

    # set length of the prefix and suffix extensions
    #if ref_species == "Hs":
    #    length = 10000   # testing from 30000 07/05/2022
    #else:
    #    length = 10000

    length = 10000
    # GUI is used to run FREEDA
    if all_genes_dict:
        # check if long exons are expected (overrides default length)
        if all_genes_dict[gene][1] is True:
            length = 30000

    prefix = 0
    longest_prefix = prefix
    while prefix < start and prefix < length:
        prefix += 1
        longest_prefix = prefix
    suffix = 0
    longest_suffix = suffix
    while suffix < seq_length and suffix < length:
        suffix += 1
        longest_suffix = suffix

    return longest_prefix, longest_suffix


def generate_files_to_MSA(contig, cds, gene, fasta_path):
    """Generates the file to align, outputs the path to that file"""

    MSA_path = fasta_path + "/MSA/"

    # select sequences to assemble an MSA file
    sequences, comp = select_contigs_to_MSA(contig, fasta_path)
    # name the MSA file depending on if it carries "comp" contig or not
    if comp is False:
        name = "to_align_" + str(contig) + ".fasta"
    else:
        name = "to_align_rev_comp_" + str(contig) + ".fasta"
    with open(name, "w") as o:
        o.write(cds.rstrip("\n"))
        o.write("\n" + sequences[0].rstrip("\n"))
        o.write("\n" + sequences[1].rstrip("\n"))
        o.write("\n" + gene.rstrip("\n"))

    o.close()
    shutil.move(name, MSA_path)

    return MSA_path


def select_contigs_to_MSA(contig, fasta_path):
    """Selects contig sequence for multiple sequence alignment from fasta folder above threshold"""

    selected_seqs = []
    rev_comp_seqs = []
    # marker if match was on the opposite strand
    comp = False

    # get paths to fasta files for a given contig
    pattern = "*_" + str(contig) + ".fasta"
    sorted_paths = sorted(glob.glob(fasta_path + "/" + pattern), key=os.path.getsize,
                          reverse=False)  # ADDED sorting (07/06/2021)

    for path in sorted_paths:

        # if match on opposite strand, fil rev_comp_seqs list with comp sequences
        if "comp" in path:
            comp = True
            with open(path, "r") as o:
                seq = o.read()
                rev_comp_seqs.append(seq)

        # if a given contig wasnt matched on the opposite strand, fill the other list
        elif comp is False:
            with open(path, "r") as o:
                seq = o.read()
                selected_seqs.append(seq)
        # if a given path contains contig that has been reversed, do nothing
        else:
            pass

    # return "comp" contig if it has been reversed, do nothing with the rest
    if comp is True:
        return rev_comp_seqs, comp

    # if contig wasnt reversed, return it
    else:
        return selected_seqs, comp
