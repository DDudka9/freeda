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

Runs MAFFT using a BioPython wrapper.

"""

from freeda import pyinstaller_compatibility
from Bio.Align.Applications import MafftCommandline
from Bio.Application import ApplicationError
import os
import re
import shutil
import glob
import time
import logging

logging.basicConfig(level=logging.INFO, format="%(message)s")


def run_msa(MSA_path):
    """Runs a multiple sequence alignment of ref cds, ref gene, presumptive locus and raw blast matches.
    Uses MAFFT as default."""

    # get path to all separate MSA files
    for in_filename in glob.glob(MSA_path + "to_align*.fasta"):

        # check if its a rev_comp file
        if re.search(r"to_align_rev_comp", in_filename):
            # for each MSA path find contig name; make it a string with group method
            out_filename = "aligned_rev_comp_" + re.search(r"(?<=to_align_rev_comp_).*$", in_filename).group()

        elif re.search(r"to_align_", in_filename):
            # for each MSA path find contig name; make it a string with group method
            out_filename = "aligned_" + re.search(r"(?<=to_align_).*$", in_filename).group()

        # run msa and record standard output and standard error
        start_time = time.time()
        handle = in_filename.split("/")[-1].replace("to_align_", "").replace("rev_comp_", "").replace(".fasta", "")
        message = "\n(%s) Aligning contig : %s with cds and gene from reference species... " \
                  % ("MAFFT", handle)
        print(message)
        logging.info(message)

        try:
            cline = MafftCommandline(cmd=pyinstaller_compatibility.resource_path("mafft"),
                                         input=in_filename,
                                         thread=-1)  # thread -1 is suppose to automatically calculate physical cores
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

        except ApplicationError as e:
            print(e)
            logging.info(e)
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
