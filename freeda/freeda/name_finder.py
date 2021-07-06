#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:25:47 2021

@author: damian

Uses regular expression module to extract protein and genome names from blast output

"""

import re

def get_names(match_path):
    
    # isolate flag names for protein and genome from blast result filename:
    # get blast result name
    path_split = re.split(r"/", match_path)[-1]
    # get protein name
    protein_name = re.split(r"_", path_split)[0]
    # get genome name (in 3 steps)
    genome_file_name = re.split(r"_", path_split)[1:3]
    genome_suffix = re.split(r"\.", genome_file_name[1])[0]
    genome_name = genome_file_name[0] + "_" + genome_suffix
    return protein_name, genome_name