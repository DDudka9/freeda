#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 18:09:48 2021

@author: damian

Runs MAFFT using a BioPython wrapper.

"""

from freeda import fasta_reader
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import MSAProbsCommandline
from Bio.Align.Applications import ProbconsCommandline
from Bio.Align.Applications import PrankCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.Application import ApplicationError
import os
import re
import shutil
import glob
import time
import logging
import subprocess


def run_msa(MSA_path, aligner):
    """Runs a multiple sequence alignment of ref cds, ref gene, presumptive locus and raw blast matches.
    Uses MAFFT as default."""

    # get path to all separate MSA files 
    for in_filename in glob.glob(MSA_path + "to_align*.fasta"):

        #if aligner == "clustalw":
        #    cline = ClustalwCommandline("clustalw2", infile=in_filename)

        # check if its a rev_comp file
        if re.search(r"to_align_rev_comp", in_filename):
            # for each MSA path find contig name; make it a string with group method
            out_filename = "aligned_rev_comp_" + re.search(r"(?<=to_align_rev_comp_).*$", in_filename).group()

        elif re.search(r"to_align_", in_filename):
            # for each MSA path find contig name; make it a string with group method
            out_filename = "aligned_" + re.search(r"(?<=to_align_).*$", in_filename).group()

        #mafft_cline = MafftCommandline(input=in_filename)

        # run msa and record standard output and standard error
        start_time = time.time()
        handle = in_filename.split("/")[-1].replace("to_align_", "").replace("rev_comp_", "").replace(".fasta", "")
        message = "\n(%s) Aligning contig : %s with cds and gene from reference species... " \
                  % (aligner.upper(), handle)
        print(message)
        logging.info(message)

        try:

            # define which aligner is used
            if aligner == "mafft":

                cline = MafftCommandline(input=in_filename, thread=-1)  # thread -1 is suppose to automatically
                                                                        # calculate physical cores
                stdout, stderr = cline()
                stop_time = time.time()
                message = "Done : in %s minutes" % ((stop_time - start_time) / 60)
                print(message)
                logging.info(message)

                # make a post-MSA file using out_filename
                with open(out_filename, "w") as f:
                    f.write(stdout)

                # move the file to MSA_path
                shutil.move(out_filename, MSA_path)

            if aligner == "muscle":  # due to crashing muscle runs at only 1 iteration !

                cmd = ['muscle', "-in", in_filename, "-quiet", "-maxiters", "2", "-out", MSA_path + out_filename]
                subprocess.call(cmd)

                # need to reorder seqs post msa
                try:
                    fasta_reader.reorder_alignment(in_filename, MSA_path + out_filename)

                except FileNotFoundError:
                    message = "...WARNING... : Aligner failed at file %s (probably too big)" % in_filename
                    print(message)
                    logging.info(message)

                stop_time = time.time()
                message = "Done : in %s minutes" % ((stop_time - start_time) / 60)
                print(message)
                logging.info(message)

                return

        except ApplicationError:
            message = "...WARNING... : Aligner failed at file %s (probably too big)" % in_filename
            print(message)
            logging.info(message)

            # rename the unaligned file as "failed"
            os.rename(in_filename, in_filename.replace("to_align", "FAILED_to_align"))

            stop_time = time.time()
            message = "FAILED : in %s minutes" % ((stop_time - start_time) / 60)
            print(message)
            logging.info(message)

            return

