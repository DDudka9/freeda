#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 18:09:48 2021

@author: damian

Runs MAFFT using a BioPython wrapper.

"""

from Bio.Align.Applications import MafftCommandline
import re
import shutil
import glob


def run_MAFFT(MSA_path, mafft_path):
    # get path to all separate MSA files 
    for in_filename in glob.glob(MSA_path + "to_align*.fasta"):
        # check if its a rev_comp file
        if re.search(r"to_align_rev_comp", in_filename):
            # for each MSA path find contig name; make it a string with group method
            out_filename = "aligned_rev_comp_" + re.search(r"(?<=to_align_rev_comp_).*$", in_filename).group()
        elif re.search(r"to_align_", in_filename):
            # for each MSA path find contig name; make it a string with group method
            out_filename = "aligned_" + re.search(r"(?<=to_align_).*$", in_filename).group()        
        # run mafft
        mafft_cline = MafftCommandline(mafft_path, input=in_filename)
        #print(mafft_cline)
        # record standard output and standard error
        stdout, stderr = mafft_cline()
        # make a post-MSA file using out_filename
        with open(out_filename, "w") as f:
            f.write(stdout)
            # move the file to MSA_path
            shutil.move(out_filename, MSA_path)