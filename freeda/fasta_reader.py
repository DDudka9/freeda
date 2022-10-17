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

Contains functions helpful to parse fasta files.

"""

from freeda import input_extractor
import re
import logging

logging.basicConfig(level=logging.INFO, format="%(message)s")


def reorder_alignment(in_filename, out_filename):
    """Reorders files post alignment based in the order before alignment"""

    ord_headers = []

    # read sequence order before the alignment
    with open(in_filename, "r") as f:
        file = f.readlines()
        for line in file:
            if line.startswith(">"):
                ord_headers.append(line.rstrip("\n"))

    headers = []
    seqs = []
    seq = ""
    count = 0

    # read the sequence order after alignment
    with open(out_filename, "r") as f:
        file = f.readlines()
        for line in file:
            count += 1

            if line.startswith(">") and count == 1:
                headers.append(line.rstrip("\n"))

            elif not line.startswith(">") and len(headers) == 1:
                seq += line.rstrip("\n")

            elif line.startswith(">") and count > 1:
                seqs.append(seq)
                headers.append(line.rstrip("\n"))
                seq = ""

            else:
                seq += line.rstrip("\n")

        # append last sequence
        seqs.append(seq)

    aln_dict = {}
    for head in headers:
        aln_dict[head] = seqs.pop(0)

    ord_aln_dict = {}
    for head in ord_headers:
        ord_aln_dict[head] = aln_dict[head]

    # overwrite the msa file with one of correct order
    with open(out_filename, "w") as f:
        for head, seq in ord_aln_dict.items():
            f.write(head + "\n")
            f.write(seq + "\n")


def alignment_file_to_dict(wdir, ref_species, filename, raw_cds=False):
    """Reads alignment file and transforms it into python dictionary.
    Alignment needs to be in the working directory (Data folder)."""

    # need to cancel wdir for path compliance
    if raw_cds:
        wdir = ""

    all_seq_dict = {}
    seq = ""

    # read alignment file and put into dict
    with open(wdir + filename, "r") as f:
        file = f.readlines()

        for line in file:

            if ">" in line and ref_species in line:
                header = line.rstrip("\n")
                continue

            if not line.startswith(">"):
                seq += line.rstrip("\n")
                continue

            # adds ref species header the first time its executed
            if line.startswith(">"):
                all_seq_dict[header] = seq
                header = line.rstrip("\n")
                seq = ""
                continue

        # record the final species
        all_seq_dict[header] = seq

    return all_seq_dict


def read_fasta_record(record):
    """Splits a fasta file into a header and sequence."""

    header = ">"
    seq = ""
    for line in re.split("\n", record):
        if line.startswith(">"):
            head = line.lstrip(">").rstrip("\n")
            header = header + "_" + head
        else:
            seq = seq + line.rstrip("\n").upper()
    return header, seq


def find_gene_and_cds(wdir, gene, ref_species):
    """Reads exons, CDS and gene of reference species making sure CDS is a one-liner"""

    ref_exons, expected_exons = get_ref_exons(wdir, gene, ref_species)

    # open according cds fasta file
    with open(wdir + "Coding_sequences/" + gene + "_" + ref_species + "_cds.fasta", "r") as f:
        make_linear = [line.rstrip("\n") for line in f.readlines()]
        for line in make_linear:
            if line.startswith(">"):
                cds_linear = line + "\n"
            else:
                cds_linear += line

    cds = cds_linear

    # open according gene fasta file
    with open(wdir + "Genes/" + gene + "_" + ref_species + "_gene.fasta", "r") as f:
        gene = f.read()

    return cds, gene, ref_exons, expected_exons


def get_ref_exons(wdir, gene, ref_species, at_input=False):
    """Reads reference species exons from input exons file into a dict used to cloned each single exon from MSA"""

    # get path to the exons for given gene
    ref_exons = {}
    seq_recorded = False
    header = ""
    seq = ""

    with open(wdir + "Exons/" + gene + "_" + ref_species + "_exons.fasta", "r") as f:
        file = f.read()

        for line in re.split("\n", file):

            # this statement executes last
            if line.startswith(">") and seq_recorded is True:
                # record header, sequence and length of the exon
                ref_exons[nr] = (header, seq, len(seq))
                seq_recorded = False
                seq = ""

            # this statement executes first
            if line.startswith(">") and seq_recorded is False:
                nr = int(line.split("_")[-1])
                head = line.lstrip(">").rstrip("\n")
                header = ">" + head
                seq_recorded = True

            # this statement executes next
            if not line.startswith(">"):
                seq = seq + line.rstrip("\n").upper()

        # record the last exon
        ref_exons[nr] = (header, seq, len(seq))

        expected_exons = tuple(e for e, features in ref_exons.items())

        # dont print and log it if function used by input extractor module
        if at_input is False:
            message = "........................................................................\n\n" \
                      "ANALYZING GENE: %s \n\n" \
                      "........................................................................\n\n" \
                      "Expected exons : %s" % (gene, str(expected_exons))
            logging.info(message)

        # flag potential microexons
        microexons = input_extractor.check_microexons(wdir, gene, ref_species)

        # dont print and log it if function used by input extractor module
        if at_input is False:
            if microexons:
                message = "\n...WARNING... : Exon nr %s is a microexon -> hard to align -> eliminated\n" % microexons
                logging.info(message)

    return ref_exons, expected_exons