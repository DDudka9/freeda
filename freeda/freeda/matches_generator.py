#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:32:54 2021

@author: damian

Uses pandas to parse blast output into a dataframe with location
of genomic loci of interest

"""

import pandas as pd
import logging

def generate_matches(match_path, t, protein_name, genome_name, genome_index):
    columns = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send",
               "evalue", "bitscore", "length", "pident", "mismatch", 
               "gapopen", "qlen", "slen"]
    
    # parse this blast table into a dataframe
    m_dataframe = pd.read_csv(match_path, sep="\t", names=columns)
    # assign variable types ("object" is default in dataframe)
    m_assigned_dataframe = m_dataframe.astype({"sseqid": str,
                        "sstart": int,
                        "send": int,
                        "pident": float})
    # select matches above treshold (30% defult) and add template "strand" Series
    selected_matches = threshold_matches(m_assigned_dataframe, t, protein_name, genome_name)
    # sort matches by "sstart" column to be able to split indicated contigs
    sorted_selected_matches = selected_matches.sort_values(by=["sseqid","sstart"])
    # reindex the rows to itertare over them easier
    matches = sorted_selected_matches.reset_index(drop=True)
    
    # record all contigs and starts/ends
    d = make_contigs_dict(matches)
    dataframes = split_contigs_to_dataframes(matches, d)
    # split matches if further than most known mouse genes (>30kb)
    concatenated_matches = split_large_contigs(dataframes).reset_index(drop=True)
    matches = concatenated_matches.astype({"sstart": int, "send": int})
    
    # record all contigs and starts/ends for strand determination
    d2 = make_contigs_dict(matches)
    dataframes2 = split_contigs_to_dataframes(matches, d2)
    # determine on which strand a given match is and split to diffrent contigs
    concatenated_matches = assign_strand(dataframes2).reset_index(drop=True)
    matches = concatenated_matches.astype({"sstart": int, "send": int})
    
    # record all contigs and starts/ends for strand splitting
    d3 = make_contigs_dict(matches)
    dataframes3 = split_contigs_to_dataframes(matches, d3)
    # split opposite strands to new contigs
    concatenated_matches = split_strands(dataframes3).reset_index(drop=True)
    matches = concatenated_matches.astype({"sstart": int, "send": int})
    
    # sort by contig name and start
    sorted_matches = concatenated_matches.sort_values(by=["sseqid","sstart"])
    matches = sorted_matches.astype({"sstart": int, "send": int}).reset_index(drop=True)
    
    return matches


def threshold_matches(matches, t, protein_name, genome_name): # works well
    # add empty "strand" column
    strand = ["for" for index in range(len(matches))]
    matches["strand"] = strand
    # filter matches based on identity 
    threshold = matches["pident"] > t
    
    # get only columns of interest and check that there is less than 40 matches (ease on MAFFT)
    matches_above_threshold = matches[threshold]
    if len(matches_above_threshold) < 40:
        pass
    
    # if not, increase identity threshold until getting < 40 matches
    else:
        while len(matches_above_threshold) > 40:
            t = t + 5
            threshold = matches["pident"] > t
            m = matches[threshold]
            matches_above_threshold = m

    side_note = "-------------------------------------------------\n" \
        "\n----- Protein: " + protein_name + " in genome: " + genome_name + \
        " -> matches identity threshold used: "+ str(t) + "\n" \
        "\n-------------------------------------------------\n"
    print(side_note)
    logging.info(side_note)
        
    # get rid of duplicated contig names
    selected_matches = matches_above_threshold[["sseqid", "sstart", \
                                                "send", "qseqid", "strand"]]
    return selected_matches


def make_contigs_dict(matches):
    d = {}
    for i, j in matches.iterrows():
        # get whole row and sstart value
        row, start, end, strand = matches.iloc[i], j[1], j[2], j[3]
        # get contig name
        contig_name = row["sseqid"]
        # get the index from Series a
        contig_idx = row.name
        # write sstart and sseqid tuple into the d dictionary under index for this row
        d[contig_idx] = contig_name, start, end, strand
    return d


def split_contigs_to_dataframes(matches, d):
    # collect contig names from matches dictionary
    contig_names = set(contig_name[0] for idx, contig_name in d.items())
    # make a list of all contig names with matches as values
    dataframes = []
    for contig_name in contig_names:
        m = matches[matches["sseqid"] == contig_name]
        m = m.reset_index(drop=True)
        dataframes.append(m) 
    return dataframes


def split_large_contigs(dataframes):
    # make a list of new matches, split into new contig names if larger than most mouse introns (30kb)
    # Long and Deutsch 1999 Nucleic Acids Research
    
    list_new_matches = []
    for matches in dataframes:
        number = 1
        # make an empty matches dataframe template to avoid copying the orginal
        new_matches = pd.DataFrame().reindex_like(matches)
        # populate the new dataframe with new contig names whenever > 30kb
        for index, row in matches.iterrows():
            start = row[1]
            if index == 0:
                contig_name = row[0]
                new_matches.iloc[index] = matches.iloc[index]
                continue
            if start - new_matches.iloc[index-1][1] < 30000:
                new_matches.iloc[index] = matches.iloc[index]
                new_matches.at[index, "sseqid"] = contig_name
                continue
            if start - new_matches.iloc[index-1][1] > 30000:
                contig_name = row[0] + "__" + str(number)
                number += 1
                new_matches.iloc[index] = matches.iloc[index]
                new_matches.at[index, "sseqid"] = contig_name
                continue
            
        list_new_matches.append(new_matches)    
    # concatenate all the new dataframes into one
    concatenated_matches = pd.concat([i for i in list_new_matches])
    return concatenated_matches


def assign_strand(dataframes2):
    # make a list of new matches, each match will have a "strand" defined \
    # forward if "0" or reverse if "1" (for a given contig)
    list_new_matches = [] 
    for matches in dataframes2:
        # make an empty matches dataframe template to avoid copying the orginal
        new_matches = pd.DataFrame().reindex_like(matches)
        for index, row in matches.iterrows():
            start, end = row[1], row[2]
            if start > end:
                new_matches.iloc[index] = matches.iloc[index]
                new_matches.at[index, "strand"] = "rev"
            else:
                new_matches.iloc[index] = matches.iloc[index]
        list_new_matches.append(new_matches)    
    # concatenate the matches into one new dataframe
    concatenated_matches = pd.concat([i for i in list_new_matches])
    return concatenated_matches


def split_strands(dataframes3):
    # make a list of new matches, split into new contig dependent on strand
    list_new_matches = [] 
    for matches in dataframes3:
        # make an empty matches dataframe template to avoid copying the orginal
        new_matches = pd.DataFrame().reindex_like(matches)
        for index, row in matches.iterrows():
            strand = row[4]
            contig_name = row[0] + "__" + strand
            new_matches.iloc[index] = matches.iloc[index]
            new_matches.at[index, "sseqid"] = contig_name
        list_new_matches.append(new_matches)    
    # concatenate the matches into one new dataframe
    concatenated_matches = pd.concat([i for i in list_new_matches])
    return concatenated_matches








