#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 19:01:01 2021

@author: damian

Splits a fasta file into a header and sequence.

"""
import re


def read_fasta_record(record):
    header = ">"
    seq = ""
    for line in re.split("\n", record):
        if line.startswith(">"):
            head = line.lstrip(">").rstrip("\n")
            header = header + "_" + head
        else:
            seq = seq + line.rstrip("\n").upper()
    return header, seq