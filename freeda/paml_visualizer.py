#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  9 14:12:04 2021

@author: damian
"""

from freeda import paml_launcher
from freeda import fasta_reader

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl import Workbook
from openpyxl.styles import Border, Side, PatternFill, Font, Alignment
from json import dump

import re
import os
import shutil
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math


wdir = "/Volumes/DamianEx_2/Data/"
result_path = wdir + "Results-08-23-2021-00-24/"
all_proteins = ["CD46"]
nr_of_species_total_dict = {"CD46" : 9}
PAML_logfile_name = "PAML-08-23-2021-00-53.log"
original_species = "Hs"
day = "-06-27-2021-15-04"
protein_name = "CD46"
nr_of_species_total = 9


# ADD LEGEND BOX TO PYMOL SCRIPT -> pymol.cgo module


def analyse_PAML_results(wdir, result_path, all_proteins, nr_of_species_total_dict,
                         original_species, PAML_logfile_name, day):
    
    all_matched_adaptive_sites_original = {}
    
    for protein in all_proteins:
        
        protein_folder = result_path + protein + "/"
        if os.path.exists(protein_folder + protein + "_PAML_analysis.tif"):
            message = "\n***************\n\n PAML analysis graph for : %s already exists" % protein
            print(message)
            logging.info(message)
        
        else:
            nr_of_species_total = nr_of_species_total_dict[protein]
            matched_adaptive_sites_original = plot_PAML(wdir, result_path, protein, nr_of_species_total, original_species)
            all_matched_adaptive_sites_original[protein] = matched_adaptive_sites_original
        
    # prepare a dict with PAML stats
    final_PAML_log_dict = read_output_PAML(result_path, PAML_logfile_name, all_matched_adaptive_sites_original)
    # generate a PAML result excel sheet
    output_excel_sheet(wdir, final_PAML_log_dict, result_path, day)
    
    # save the matched adaptive sites into a file
    with open("all_matched_adaptive_sites_original.txt", "w") as file:
        dump(all_matched_adaptive_sites_original, file)
    
    shutil.move("all_matched_adaptive_sites_original.txt", result_path)


def read_output_PAML(result_path, PAML_logfile_name, all_matched_adaptive_sites_original):

    # prepare the dict storing the PAML result
    PAML_log_dict = {"Protein name":[],
                        "Nr of species analyzed":[],
                        "Species":[],
                        "CDS Coverage":[],
                        "M2a vs M1a (LRT)":[],
                        "M2a vs M1a (p-value)":[],
                        "M8 vs M7 (LRT)":[],
                        "M8 vs M7 (p-value)":[],
                        "Sites with pr < 0.90":[],
                        "Sites with pr >= 0.90":[]} 
    
    # extract PAML results from PAML log file
    PAML_log_path = result_path + PAML_logfile_name
    
    with open(PAML_log_path, "r") as f:
        file = f.readlines()
        
        start_recording = False

        for line in file:
            
            # record protein name
            if line.startswith(" --------- *"):
                protein = line.replace("-","").replace(" ","").replace("*","").rstrip("\n")
                PAML_log_dict["Protein name"].append(protein)
                start_recording = True 

            # record species
            if start_recording is True and line.startswith(" Final species"):
                species = line.split(":")[1].split("[")[1].replace("]", "").replace("\n","").replace("'", "")
                nr_species = line.split(":")[1].split(" ")[1]
                PAML_log_dict["Nr of species analyzed"].append(nr_species)
                PAML_log_dict["Species"].append(species)
                
            # record CDS_coverage
            if start_recording is True and line.startswith("['Original alignment:"):
                coverage = line.split(")")[0].split("(")[-1].replace(" ", "")
                PAML_log_dict["CDS Coverage"].append(coverage)
            
            # record LRTs
            if start_recording is True and line.startswith(" PAML LRTs"):
                
                LRT1 = line.split(":")[1].split("and")[0].split("-")
                if LRT1[1].endswith("e") is True:
                    M2a_vs_M1a_LRT = LRT1[1].replace(" ", "") + "-" + LRT1[2]
                if LRT1[1].endswith("e") is False:
                    M2a_vs_M1a_LRT = LRT1[1].replace(" ", "")
                if M2a_vs_M1a_LRT != "None":
                    M2a_vs_M1a_LRT = round(float(M2a_vs_M1a_LRT), 4)
                    if M2a_vs_M1a_LRT == 0.0:
                        M2a_vs_M1a_LRT = 0.0001

                LRT2 = line.split(":")[1].split("and")[1].split("-")
                # catch "e-" notation because it gets split
                if LRT2[1].endswith("e") is True:
                    M8_vs_M7_LRT = LRT2[1].replace(" ", "") + "-" + LRT2[2]
                if LRT2[1].endswith("e") is False:
                    M8_vs_M7_LRT = LRT2[1].replace(" ", "")
                if M8_vs_M7_LRT != "None":
                    M8_vs_M7_LRT = round(float(M8_vs_M7_LRT), 4)
                    if M8_vs_M7_LRT == 0.0:
                        M8_vs_M7_LRT = 0.0001
                    
                PAML_log_dict["M2a vs M1a (LRT)"].append(M2a_vs_M1a_LRT)
                PAML_log_dict["M8 vs M7 (LRT)"].append(M8_vs_M7_LRT)
            
            # record p_values
            if start_recording is True and line.startswith(" -> PAML p-values"):
                
                p_value1 = line.split(":")[1].split("and")[0].split("-")
                # catch "e-" notation because it gets split
                if p_value1[1].endswith("e") is True:
                    M2a_vs_M1a_pvalue = p_value1[1].replace(" ","") + "-" + p_value1[2]
                if p_value1[1].endswith("e") is False:
                    M2a_vs_M1a_pvalue = p_value1[1].replace(" ","")
                if M2a_vs_M1a_pvalue != "None":
                    M2a_vs_M1a_pvalue = round(float(M2a_vs_M1a_pvalue), 4)
                    if M2a_vs_M1a_pvalue == 0.0:
                        M2a_vs_M1a_pvalue = 0.0001
                    
                PAML_log_dict["M2a vs M1a (p-value)"].append(M2a_vs_M1a_pvalue)
                
                p_value2 = line.split(":")[1].split("and")[1].split("-")
                # catch "e-" notation because it gets split
                if p_value2[1].endswith("e"):
                    M8_vs_M7_pvalue = p_value2[1].replace(" ","") + "-" + p_value2[2].rstrip("\n")
                if p_value2[1].endswith("e") is False:
                    M8_vs_M7_pvalue = p_value2[1].replace(" ","").rstrip("\n")
                if M8_vs_M7_pvalue != "None":
                    M8_vs_M7_pvalue = round(float(M8_vs_M7_pvalue), 4)
                    if M8_vs_M7_pvalue == 0.0:
                        M8_vs_M7_pvalue = 0.0001

                PAML_log_dict["M8 vs M7 (p-value)"].append(M8_vs_M7_pvalue)
                    
                # does not report adaptive sites in proteins failing M8 vs M7 test
                if M8_vs_M7_pvalue == "None" or float(M8_vs_M7_pvalue) >= 0.05:
                    PAML_log_dict["Sites with pr < 0.90"].append("Not likely")
                    PAML_log_dict["Sites with pr >= 0.90"].append("Not likely")
                
                # remove adaptive sites from proteins that are not rapidly evolving
                if M8_vs_M7_pvalue != "None" and float(M8_vs_M7_pvalue) <= 0.05:
                        
                    matched_dict = all_matched_adaptive_sites_original[protein]    
                    mild_sites = {}
                    strong_sites = {}
        
                    for site, values in matched_dict.items():

                        if float(values[2]) > 0.00:
                            adaptive_site = values[0] + str(site)
                            probability = float(values[2])
                                
                            if float(values[2]) < 0.90:
                                mild_sites[adaptive_site] = probability
    
                            if float(values[2]) >= 0.90:
                                strong_sites[adaptive_site] = probability
                        
                    mild_sites_to_append = (str(len(mild_sites)) + " " + str(mild_sites)).replace("'","")
                    strong_sites_to_append = (str(len(strong_sites)) + " " + str(strong_sites)).replace("'","")
                        
                    PAML_log_dict["Sites with pr < 0.90"].append(mild_sites_to_append)
                    PAML_log_dict["Sites with pr >= 0.90"].append(strong_sites_to_append)

                # finish recording loop through the log file
                start_recording = False
                
        final_PAML_log_dict = PAML_log_dict
        
    return final_PAML_log_dict


def output_excel_sheet(wdir, final_PAML_log_dict, result_path, day):

    # make an empty excel sheet using openpyxl library
    wb = Workbook()
    ws = wb.active
    
    # set column width
    columns_to_format = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    for c in columns_to_format:
        ws.column_dimensions[c].width = 19

    # format headers (cell by cell)
    headers = [ws["A1"], ws["B1"], ws["C1"], ws["D1"], ws["E1"], 
               ws["F1"], ws["G1"], ws["H1"], ws["I1"], ws["J1"]]
    heads = list(final_PAML_log_dict.keys())
    header_count = 1
    for header in headers:
        
        # THIS IS MOST LIKELY NOT NEEDED (CHECK)
        head = heads.pop(0)
        ws.cell(row=1, column=header_count, value=head)
        
        header.font = Font(bold = True, name = "Arial", size = 10)
        double = Side(border_style = "thin")
        header.border = Border(top = double, left = double, right = double, bottom = double)
        header.alignment = Alignment(horizontal = "center", vertical = "center")
        header.fill = PatternFill("solid", fgColor = "00FFCC00")    
        header_count += 1
    
    all_rows = []
    for entry in final_PAML_log_dict.values():
    
        # make it a dataframe
        df = pd.DataFrame(entry)

        # read in all values from each key of final_PAML_log_dict (remove index column)
        rows = dataframe_to_rows(df, index=False)
        all_rows.append(list(rows))
    
    # counters for the iteration over cells in the worksheet (openpyxl index starts at 1)
    # leave the first row for headers
    column_count = 1
    row_count = 2
    
    # go over each entry per list of entries from the final PAML log dict
    for row in all_rows:

        for one_entry in row:
            
            # need to convert the generator into a string to tidy it
            entry_as_string = str(one_entry).replace("[", "").replace("]", "").replace("'", "").replace("\n", "")
            one_entry = entry_as_string
            
            # reset the row counter
            if row_count == len(list(row))+1:
                column_count += 1
                row_count = 2
            
            # record a header of new column
            if one_entry == "0":
                head = [key for number, key in enumerate(final_PAML_log_dict.keys(),1) if number == column_count][0]
                
                # always goes into the first row
                ws.cell(row=1, column=column_count, value=head.replace("'",""))
            
            # record entry
            else:
                ws.cell(row=row_count, column=column_count, value=one_entry)
                row_count += 1
    
    # make all rows centered
    for row in ws.iter_rows():
        for cell in row:
            cell.alignment = Alignment(horizontal = "center", vertical = "center")
    
    # save the excel document
    excel_filename = "PAML_result" + day + ".xlsx"
    wb.save(excel_filename)
    shutil.move(wdir + excel_filename, result_path)



def plot_PAML(wdir, result_path, protein, nr_of_species_total, original_species):
    
    # get original and final aa sequence for a protein
    original_sequence_record, final_sequence_record = get_original_and_final_seqs(wdir, protein, result_path, original_species)
    # organise them into easy to search dictionaries
    original_species_dict, final_original_species_dict, final_length = organise_original_and_final_seqs(original_sequence_record, final_sequence_record)
    # map the aa residues between the original and final sequnces
    mapped_original_and_final_residues_dict = map_original_and_final_residues(original_sequence_record, final_sequence_record)
    # find dN/dS omega ratio per site based on "rst" file (final only for now)
    omega_dict = get_omegas(protein, result_path, final_length)
    # read PAML output file to find which residues are rapidly evolving in the final seq
    adaptive_sites_dict = get_adaptive_sites(result_path, protein)
    # find these residues in the original sequence
    matched_adaptive_sites_original = match_adaptive_sites_to_original(final_original_species_dict,
            mapped_original_and_final_residues_dict, adaptive_sites_dict, omega_dict)
    # mark the sites that were not analyzed by PAML
    final_dict_to_plot = mark_skipped_sites(matched_adaptive_sites_original, mapped_original_and_final_residues_dict)
    # record and write breakdown of adaptive sites overlay to original cds
    record_adaptive_sites(final_dict_to_plot, protein)
    # plot omegas and probabilities
    make_graphs(final_dict_to_plot, result_path, protein, nr_of_species_total)

    return matched_adaptive_sites_original

def mark_skipped_sites(matched_adaptive_sites_original, mapped_original_and_final_residues_dict):
    
    # mark the sites that were not analyzed by PAML    
    for site in mapped_original_and_final_residues_dict:
        if mapped_original_and_final_residues_dict[site][0] != "-":
            matched_adaptive_sites_original[site].append(1)
        else:
            matched_adaptive_sites_original[site].append(0)
            
    return matched_adaptive_sites_original


def record_adaptive_sites(final_dict_to_plot, protein_name):
    
    # row_residues and row_features are well set
    row_positions = ""
    row_residues = ""
    row_features = ""
    position_counter = 0
    row_counter = [row for row in range(1, len(final_dict_to_plot) + 1) if row % 10 == 0]
    
    for position, values in final_dict_to_plot.items():
        
        # break all rows for visibility
        if position_counter % 50 != 0 and position_counter != 0 and position_counter % 10 == 0:

            row_residues = row_residues + " "
            row_features = row_features + " "
        
        position_counter += 1
        
        # prepare a position count array
        # first row
        if position_counter == 1:
            position_marker = row_counter.pop(0)
            row_positions = 8 * " " + str(position_marker) + 9 * " "
        
        # any other row
        if (position_counter + 1) % 10 == 0 and row_counter != []:
            position_marker = row_counter.pop(0)
            
            if position_marker < 100:
                
                if position_marker % 50 == 0:
                    row_positions = row_positions + str(position_marker) + "\n" + 8 * " "
                else:
                    row_positions = row_positions + str(position_marker) + 9 * " "
                    
            if position_marker == 100:
                row_positions = row_positions[:-1] + str(position_marker) + "\n" + 7 * " "
                    
            if 1000 > position_marker > 100:
        
                if position_marker % 50 == 0:
                    row_positions = row_positions + str(position_marker) + "\n" + 7 * " "
                else:
                    row_positions = row_positions + str(position_marker) + 8 * " "
            
            if position_marker == 1000:
                row_positions = row_positions[:-1] + str(position_marker) + "\n" + 6 * " "
            
            if 1000 < position_marker:
            
                if position_marker % 50 == 0:
                    row_positions = row_positions + str(position_marker) + "\n" + 6 * " "
                else:
                    row_positions = row_positions + str(position_marker) + 7 * " "
        
        # prepare residues and features arrays
        # not adaptive
        if float(values[2]) < 0.75 and values[3] != 0:
            residue = values[0]
            row_features = row_features + " "
            row_residues = row_residues + residue
            if position % 50 == 0:
                row_residues = row_residues + "\n"
                row_features = row_features + "\n"
            
        # mild probability of adaptive evolution
        if 0.75 <= float(values[2]) < 0.90:
            residue = values[0]
            row_features = row_features + "."
            row_residues = row_residues + residue
            if position % 50 == 0:
                row_residues = row_residues + "\n"
                row_features = row_features + "\n"
                
        # strong probability of adaptive evolution
        if 0.90 <= float(values[2]):
            residue = values[0]
            row_features = row_features + ":"
            row_residues = row_residues + residue
            if position % 50 == 0:
                row_residues = row_residues + "\n"
                row_features = row_features + "\n"
        
        # residue absent in PAML analysis
        if values[3] == 0:
            residue = values[0]
            row_features = row_features + "-"
            row_residues = row_residues + residue
            if position % 50 == 0:
                row_residues = row_residues + "\n"
                row_features = row_features + "\n"

    positions = row_positions.split("\n")
    features = row_features.split("\n")
    residues = row_residues.split("\n")
    
    if len(positions) == len(features) == len(residues):
    
        message = ("\n\n.........................................."
            "\n\nReference sequence for %s with adaptive sites:"
            "\n\n . means pr >= 0.75 \n : means pr >= 0.90 \n - means missing from PAML analysis\n\n") % protein_name
        print(message)
        logging.info(message)
        
        for row in range(len(positions)):
            message = positions[row]
            print(message)
            logging.info(message)
            
            message = features[row]
            print(message)
            logging.info(message)
            
            message = residues[row]
            print(message)
            logging.info(message)
            
    else:
        message = "problem with unequal rows"
        print(message)
        print(positions)
        print(features)
        print(residues)
        logging.info(message)
        logging.info(positions)
        logging.info(features)
        logging.info(residues)
    
    

    

def get_original_and_final_seqs(wdir, protein_name, result_path, original_species):
    
    protein_path = result_path + "/" + protein_name

    # get record onject for the iriginal protein sequence
    original_sequence = paml_launcher.get_original_cds(wdir, protein_name, original_species)
    original_sequence_Seq = Seq(original_sequence)
    original_original_Seq_object = original_sequence_Seq.translate()
    original_sequence_record = SeqRecord(original_original_Seq_object)
    
    # get record onject for the final protein sequence
    pattern = "translated.fasta"
    all_files = os.listdir(protein_path)
    
    for file in all_files:
        if pattern in file and "reduced" not in file and not file.startswith("."):
            post_Gblocks_translated_file_path = protein_path + "/" + file
    
    alignment = AlignIO.read(post_Gblocks_translated_file_path, "fasta")
    all_seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]
    final_sequence = all_seqs[0][1]
    final_sequence_Seq_object = Seq(final_sequence)
    final_sequence_record = SeqRecord(final_sequence_Seq_object)
    
    return original_sequence_record, final_sequence_record
    

def organise_original_and_final_seqs(original_sequence_record, final_sequence_record):
    
    # organise residues into easy to search dict
    original_species_dict = {}    
    total_length = 0
    for index, aa in enumerate(original_sequence_record.seq, start=1):
        if aa == "*":
            break
        original_species_dict[index] = aa
        total_length += 1
    
    # organise residues into easy to search dict
    final_original_species_dict = {}
    final_length = 0
    for index, aa in enumerate(final_sequence_record.seq, start=1):
        if aa == "*":
            break
        final_original_species_dict[index] = aa
        final_length += 1

    return original_species_dict, final_original_species_dict, final_length   


def map_original_and_final_residues(original_sequence_record, final_sequence_record):
    
    # perform pairwise alignment (multiple alignments of the same sequences)
    aln = pairwise2.align.globalxx(final_sequence_record.seq, original_sequence_record.seq)
    
    # make a dict storing the alignments
    mapped_original_and_final_residues_dict = {}
    for i in range(1, len(aln[0][0])):
        mapped_original_and_final_residues_dict[i] = []

    # fill the dict with paired positions (need to subtract 1 cose aln starts at 0)
    for position in mapped_original_and_final_residues_dict:
        mapped_original_and_final_residues_dict[position] = [aln[0][0][position-1], aln[0][1][position-1]]
    
    return mapped_original_and_final_residues_dict
    
    # flag any non-synonymous substitutions that indicate frameshifts
    #frameshift_positions = {}
    #translated_frameshift = False
    #for position, pair in d.items():
    #    if pair[0] != pair[1] and pair[0] != "-" and pair[1] != "-":
    #        # sequence index starts at 1 so need to "+1" from the python indexing
    #        frameshift_positions[position + 1] = pair
    #        translated_frameshift = True
    #        print(translated_frameshift)

    

def match_adaptive_sites_to_original(final_original_species_dict, mapped_original_and_final_residues_dict, adaptive_sites_dict, omega_dict):
    
    # overlay the final residue dictionary with probability for positive selection
    matched_adaptive_sites_final = {}
    for site, aa in final_original_species_dict.items():
        # mark this site as potentially positively selected
        if site in adaptive_sites_dict:
            matched_adaptive_sites_final[site] = adaptive_sites_dict[site]
        # mark this site as not positively selected
        else:
            matched_adaptive_sites_final[site] = [aa, "0.00"]
    
    # adjust the final residue dicitionary to include missing sites
    adjusted_mapped_original_and_final_residue_dict = {}
    for i in range(1, len(mapped_original_and_final_residues_dict) + 1):
        adjusted_mapped_original_and_final_residue_dict[i] = ["", "", 0]
    
    a = mapped_original_and_final_residues_dict
    b = adjusted_mapped_original_and_final_residue_dict
    for site, aa in a.items():
        # special case when there is no previous site/key
        if site == 1 and aa[0] == "-":
            b[site] = [aa[0], aa[1], 1]
            continue
        elif site == 1 and aa[0] != "-":
            b[site] = [aa[0], aa[1], 0]
        # collect dashes and add them to tally per site
        elif site != 1 and aa[0] == "-":
            b[site] = [aa[0], aa[1], b[site-1][2]+1]
        # propagate the tally from previous site
        elif site != 1 and aa[0] != "-":
            b[site] = [aa[0], aa[1], b[site-1][2]]
        # safety print -> REMOVE AFTER TESTING
        else:
            print("else")
    
    # map omegas and probabilities to the original residues
    c = matched_adaptive_sites_final
    w = omega_dict
    matched_adaptive_sites_original = {}
    for site, aa in a.items():
        if aa[0] == aa[1]:
            key = [k for k,v in b.items() if k == site][0]
            probability = c[key-b[site][2]][1]
            omega = w[key-b[site][2]]
            matched_adaptive_sites_original[site] = [aa[0], omega, probability]
        else:
            matched_adaptive_sites_original[site] = [aa[1], "0.000", "0.000"]

    return matched_adaptive_sites_original


def get_omegas(protein_name, result_path, final_length):
    
    PAML_output_file_path = result_path + "/" + protein_name + "/PAML_" + protein_name
    
    PAML_result_dict = {}
    for i in range(1, final_length+1):
        PAML_result_dict[i] = ""
    
    # collect all values for each residue from the BEB of M8 model
    with open(PAML_output_file_path + "/rst", "r") as f:
        
        start_search = False
        start_recording = False
        line_number = 0
        # assumes that max protein length == 9999aa
        space = " " * (4 - len(str(final_length)))
        file = f.readlines()
        
        l = 0
        for line in file:
            l += 1
            
            if line_number == final_length:
                break
            
            if start_recording is False and line.startswith("Model 8"):
                start_search = True
                continue
            
            if start_search is True and line.startswith("Bayes"):
                start_recording = True
                continue
            
            if start_recording is True and line.startswith("   1 "):
                line_number += 1
                PAML_result_dict[line_number] = line.split(" ")
                continue
            
            if start_recording is True and PAML_result_dict[1] != "":
                line_number += 1
                PAML_result_dict[line_number] = line.split(" ")
                continue
            
            if start_recording is True and line.startswith(space + str(final_length)):
                line_number += 1
                PAML_result_dict[line_number] = line.split(" ")
                start_recording = False
    
    # collect omegas for each residue
    omega_dict = {}
    for residue, values in PAML_result_dict.items():
        if float(values[-1]) > 0.001:
            omega_dict[residue] = float(values[-4])
        else:
            omega_dict[residue] = 0.000
        
    return omega_dict


def get_adaptive_sites(result_path, protein):
    
    PAML_output_file_path = result_path + "/" + protein + "/PAML_" + protein
    with open(PAML_output_file_path + "/output_PAML", "r") as f:
        
        adaptive_sites_dict = {}
        start_search1 = False
        start_search2 = False
        start_recording = False
        newlines = 0
        file = f.readlines()
        
        for line in file:

            if start_recording is True and re.search(r"^\s*\d", line):
                site_number = int(re.match(r"^\s*\d{1,4}", line).group(0).replace(" ",""))
                residue = re.match(r"^\s*\d{1,4}\s[A-Z]", line).group(0)[-1]
                probability = re.match(r"^\s*\d{1,4}\s[A-Z]\s*\d\.\d{3}", line).group(0)[-5:]
                adaptive_sites_dict[site_number] = [residue, probability]
                continue
                
            if start_recording is True and line.startswith("\n"):
                newlines += 1
                
            if newlines == 2:
                break

            if start_recording is False and line.startswith("Model 8"):
                start_search1 = True
                continue
            
            if start_search1 is True and line.startswith("Bayes"):
                start_search2 = True
                continue
            
            # 12 spaces in front of "Pr"
            if start_search2 is True and line.startswith("            Pr"):
                start_recording = True
                continue
        
    return adaptive_sites_dict


def make_graphs(final_dict_to_plot, result_path, protein, nr_of_species_total):

    sites = []
    residues = []
    omegas = []
    probabilities = []
    analyzed = []
    # these will determine omega graph y axis
    roof = 1
    floor = 1

    for site, features in final_dict_to_plot.items():
        sites.append(site)
        residues.append(features[0])
        
        if features[3] == 1:
            omegas.append(float(features[1]))
            if features[1] > roof:
                roof = math.ceil(features[1])
            elif 0.5 <= features[1] <= floor:
                floor = 0.5
            elif 0 < features[1] < 0.5:
                floor = 0
            probabilities.append(float(features[2]))
            
        # mark missing values as impossibly high omegas and probabilities
        else:
            omegas.append(float(10.00))
            probabilities.append(float(10.00))
    
        analyzed.append(features[3])
    
    # This is pretty ok
    plt.figure()

    # This is pretty ok
    plt.subplot(311, title="PAML analysis - %s (%s species analyzed)" % (protein, nr_of_species_total))
    plt.ylabel("dN/dS\n")
    # I should fix max omega better
    plt.axis([0.5, sites[-1], 0, roof])
    plt.ylim(0.5, roof)
    plt.yticks(np.arange(0.5, roof + 0.1, 1.0))
    # mark the missing values for the plot
    clrs1 = ["black" if s == 1 else "gainsboro" for s in analyzed]
    plt.bar(sites, omegas, color=clrs1)
    #plt.bar(sites, omegas, color=(0.1, 0.1, 0.1))

    
    # This is pretty ok
    plt.subplot(312)
    plt.ylabel("Prob. positive\n selection")
    plt.axis([0, sites[-1], 0.5, 1.0])
    # mark the missing values for the plot
    clrs2 = ["cornflowerblue" if s == 1 else "gainsboro" for s in analyzed]
    plt.bar(sites, probabilities, color=clrs2)
    
    # This is pretty ok
    plt.subplot(313)
    plt.xlabel("Codons (mapped to CDS used as reference)")
    plt.ylabel("High prob. \n positive selection")
    plt.axis([0, sites[-1], 0.9, 1.0])
    plt.yticks(np.arange(0.9, 1.01, 0.05))
    # mark the missing values for the plot
    clrs3 = ["magenta" if s == 1 else "gainsboro" for s in analyzed]
    plt.bar(sites, probabilities, color=clrs3)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    
    figure_name = protein + "_PAML_analysis"
    protein_path = result_path + "/" + protein
    plt.savefig(figure_name + ".tif", dpi=300, bbox_inches="tight")
    plt.savefig(figure_name + ".svg", dpi=300, bbox_inches="tight")
    
    shutil.move(figure_name + ".tif", result_path)
    shutil.move(figure_name + ".svg", protein_path)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
