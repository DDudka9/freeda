#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 20:35:28 2021

@author: damian

Generates a PyMOL script and runs it internally. Saves the PyMOL session per gene.

"""

from ast import literal_eval
import shutil
import subprocess
import os
import requests
import simplejson.errors
import sys

# LAST RESIDUE IS NOT MARKED IN PYMOL MODEL IF SCORING (C-term label interferes?)
# Done but NOT TESTED YET

# 2021-09-15
# Brian Akins
# METHOD: get_interpro(uniprot_id), interacts with the InterPro REST API.
# INPUT:  A string containing the Uniprot ID of a target gene.
# OUTPUT: A requests Response object. To interact with this object, use:
#         .status_code, .headers['content-type'], .encoding, .text, or .json()
# More info on the API URL architecture can be found at
# https://docs.google.com/document/d/1JkZAkGI6KjZdqwJFXYlTFPna82p68vom_CojYYaTAR0/edit
def get_interpro(uniprot_id):
    """Gets url from Interpro API -> json file"""
    print("\n Retrieving InterPro data for UniProt ID " + uniprot_id + "...\n")
    interpro_url = "https://www.ebi.ac.uk/interpro/api/protein/uniprot/" + uniprot_id + "/entry/interpro/"
    print("Request URL: " + interpro_url)
    response = requests.get(interpro_url)
    return response


def get_domain_info(interpro_entry_id):
    """Gets Interpro domain entry from json file"""

    domain_entry_url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/" + interpro_entry_id
    response = requests.get(domain_entry_url)
    return response


# METHOD: protein_domains(uniprot_id), gets information from the InterPro REST API about the protein domains present.
# INPUT: A string containing the UniProt ID of the target protein.
# OUTPUT: A list of dictionaries, each containing:
#         {'accession': InterPro Accession ID string,
#          'name': domain name string,
#          'coordinates': [(domain_start1, domain_end1), (domain_start2, domain_end2), ...]}
#         With one list item per InterPro entry. Most (all?) domains will have only one start/stop coordinate tuple.
def protein_domains(uniprot_id):
    """Makes python dictionary with Interpro domain entries from json file"""

    interpro_dict = get_interpro(uniprot_id).json()
    # List the entries from InterPro with type 'domain' (as opposed to 'family', etc.) as a list of dictionaries
    domains = [entry for entry in interpro_dict["entry_subset"] if entry["entry_type"] == "domain"]
    output_dict = []
    # Iterate through the domain entries and add an entry to the output dictionary for each
    for entry in domains:
        accession_id = entry["accession"]
        domain_info = get_domain_info(accession_id).json()
        domain_name = domain_info["metadata"]["name"]["name"]
        # An ugly way to retrieve the start(s) and end(s) of the domain from the nested dictionaries and lists
        # If there is a bug, check here first
        domain_fragments = entry["entry_protein_locations"][0]["fragments"]
        # A list of tuples where each tuple is (start residue number, end residue number) for the domain
        domain_coordinates = [(fragment["start"], fragment["end"]) for fragment in domain_fragments]
        output_dict.append({"accession": accession_id,
                            "name": domain_name,
                            "coordinates": domain_coordinates})
    print("Domains retrieved for UniProt ID " + uniprot_id)
    return output_dict


def get_domains(wdir, ref_species, gene):
    """Based on uniprot id gets domain layout from Interpro API : name and coordinates"""

    structure_filepath = wdir + "Structures/" + gene + "_" + ref_species

    # if input matched AlphaFold model then Uniprot ID is found in that file
    try:
        # get uniprot ID
        with open(structure_filepath + "/model_matches_input_seq.txt", "r") as f:
            file = f.readlines()
            for line in file:
                if "Uniprot ID" in line:
                    uniprot_id = line.rstrip("\n").split(":")[1].replace(" ", "")

    except FileNotFoundError:
        return

    # connection via requests seem unreliable sometimes -> 20 chances
    domains = {}
    count = 0
    while count < 20:
        try:
            interpro_dict = protein_domains(uniprot_id)
            for domain in interpro_dict:
                name = domain["name"]
                print(f"name : {name}")
                coordinates = domain["coordinates"][0]
                print(f"coordinates : {coordinates}")
                domains[name] = coordinates
            count = 20

        # no Interpro data for this gene
        except simplejson.errors.JSONDecodeError:
            return domains
        # requested url wasnt pulled correctly
        except KeyError:
            count += 1

    return domains


def check_structure(wdir, ref_species, gene):
    """Checks presence of a structure prediction model for a given protein"""

    structure_model_path = wdir + "Structures/" + gene + "_" + ref_species
    model_file_list = os.listdir(structure_model_path)
    unwanted_files = [file for file in model_file_list if file.startswith(".") or file.endswith(".png")]

    # check for unwanted files (hidden and png files) and remove them
    if unwanted_files:
        for file in unwanted_files:
            try:
                os.remove(structure_model_path + "/" + file)
            except FileNotFoundError:
                pass

    # regenerate list of files -> should contain one file exactly
    model_file_list = os.listdir(structure_model_path)

    # get info on model quality
    txt_files = [file for file in model_file_list if file.endswith(".txt")]

    # does not allow structure overlay onto incompatible model (model doesnt match input seq)
    if "model_incompatible.txt" in txt_files:
        return False

    # model matches input seq (model_incompatible.txt file takes priority)
    if "model_matches_input_seq.txt" in txt_files:
        return True

    # allow only one pdb file present
    pdb_files = [file for file in model_file_list if file.startswith(".") or not file.endswith(".pdb")]
    if len(pdb_files) == 1:
        return True

    else:
        return False


def run_pymol(wdir, ref_species, result_path, gene, genes_under_pos_sel, all_genes_dict=None):
    """Runs PyMOL with overlaid adaptive sites"""

    structures_path = wdir + "Structures"
    protein_structure_path = structures_path + "/" + gene + "_" + ref_species + "/"
    final_model_name = gene + "_" + ref_species + ".pse"

    # get dict of adaptives sites matched to input sequence
    with open(result_path + "all_matched_adaptive_sites_ref.txt", "r") as f:
        all_matched_adaptive_sites_ref = literal_eval(f.read())

    # get domain names and coordinates from Interpro API
    domains = get_domains(wdir, ref_species, gene)

    # obtain a pymol script based on the model
    if not get_pymol_script(wdir, ref_species, all_matched_adaptive_sites_ref, gene,
                            protein_structure_path, genes_under_pos_sel, domains, all_genes_dict):
        return False

    # run that script in pymol without triggering external GUI (-cq) -> DOES NOT WORK IN PYCHARM?
    pymol_command = ["pymol", "-cq", "structure_overlay.pml"]
    stderr, stdout = subprocess.Popen(pymol_command, shell=True, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).communicate()
    # move and overwrite if "structure_overlay.pml" exists in Structure folder for the gene
    shutil.move(os.path.join(wdir, "structure_overlay.pml"),
                os.path.join(protein_structure_path, "structure_overlay.pml"))
    # move the model with overlaid residues into Results folder
    shutil.move(protein_structure_path + final_model_name,
                result_path.replace("Raw_data/", "Results/Structures/") + final_model_name)

    return True


def get_consensus_dict(all_matched_adaptive_sites_ref, gene, genes_under_pos_sel):
    """Builds a consensus dictionary by comparing sites under positive selection in both codon frequency models used
    - converts sites < 0.90 into 0.00"""

    consensus_dict = {}
    # check if consensus is to be built (more than one codon frequency used)
    if len(all_matched_adaptive_sites_ref.keys()) > 1 and gene in genes_under_pos_sel["F3X4"] \
                                                        and gene in genes_under_pos_sel["F61"]:
        consensus_dict[gene] = {}
        # compare probabilities at each position between models
        for position in range(1,
                              len(all_matched_adaptive_sites_ref[next(iter(all_matched_adaptive_sites_ref))][gene]) + 1):
            position_F3X4 = all_matched_adaptive_sites_ref["F3X4"][gene][str(position)]
            position_F61 = all_matched_adaptive_sites_ref["F61"][gene][str(position)]
            # suppress scoring of residues with high probability in only one codon frequency model
            if float(position_F3X4[2]) < 0.90 or float(position_F61[2]) < 0.90:
                position_F3X4[2] = "0.00"
            consensus_dict[gene][str(position)] = position_F3X4

    # gene scores only in F3X4 model
    elif len(all_matched_adaptive_sites_ref.keys()) > 1 \
        and gene in genes_under_pos_sel["F3X4"] and gene not in genes_under_pos_sel["F61"]:

        consensus_dict[gene] = {}
        # make consensus depend only on F3X4 model
        for position in range(1,
                              len(all_matched_adaptive_sites_ref[next(iter(all_matched_adaptive_sites_ref))][
                                      gene]) + 1):
            position_F3X4 = all_matched_adaptive_sites_ref["F3X4"][gene][str(position)]
            consensus_dict[gene][str(position)] = position_F3X4

    # gene scores only in F61 model
    elif len(all_matched_adaptive_sites_ref.keys()) > 1 \
        and gene not in genes_under_pos_sel["F3X4"] and gene in genes_under_pos_sel["F61"]:

        consensus_dict[gene] = {}
        # make consensus depend only on F61 model
        for position in range(1,
                              len(all_matched_adaptive_sites_ref[next(iter(all_matched_adaptive_sites_ref))][
                                      gene]) + 1):
            position_F61 = all_matched_adaptive_sites_ref["F61"][gene][str(position)]
            consensus_dict[gene][str(position)] = position_F61

    elif len(all_matched_adaptive_sites_ref.keys()) == 1:

        consensus_dict[gene] = {}
        # only one model
        for position in range(1,
                              len(all_matched_adaptive_sites_ref[next(iter(all_matched_adaptive_sites_ref))][
                                      gene]) + 1):
            # get the first (and only) key
            positions = all_matched_adaptive_sites_ref[next(iter(all_matched_adaptive_sites_ref))][gene][str(position)]
            consensus_dict[gene][str(position)] = positions

    else:
        print("Something went wrong with the consensus dict")

    return consensus_dict


def get_pymol_script(wdir, ref_species, all_matched_adaptive_sites_ref, gene,
                     protein_structure_path, genes_under_pos_sel, domains, all_genes_dict=None):
    """Gets a PyMOL script that will be passed into PyMOL automatically"""

    paint_sites = False

    # paint sites of genes likely under positive selection
    if gene in genes_under_pos_sel["F3X4"] or gene in genes_under_pos_sel["F61"]:
        paint_sites = True

    # get consensus (more than one codon frequency used)
    consensus_dict = get_consensus_dict(all_matched_adaptive_sites_ref, gene, genes_under_pos_sel)
    matched_adaptive_sites_ref = consensus_dict[gene]

    #if consensus_dict:
    #    matched_adaptive_sites_ref = consensus_dict[gene]
    #else:
    #    matched_adaptive_sites_ref = all_matched_adaptive_sites_ref[next(iter(all_matched_adaptive_sites_ref))][gene]
                                                                                # pick the first key
                                                                                # (only one codon model)

    structure_prediction_path = wdir + "Structures/" + gene + "_" + ref_species

    if len(os.listdir(structure_prediction_path)) == 0:
        print("\nNo structure prediction model is present for: %s -> cannot run PyMOL" % gene)
        return False

    with open("structure_overlay.pml", "w") as f:

        # PyMOL command to load model
        structure_model_path = wdir + "Structures/" + gene + "_" + ref_species

        # select the model (check_structure function removed all hidden files)
        model = [file for file in os.listdir(structure_model_path) if file.endswith(".pdb")][0]

        # start pymol script by loading the model
        f.write("load " + protein_structure_path + model + "\n")

        # PyMOL command to color all residues
        f.write("color cyan\n")

        # color domains
        colors = ["orange", "marine", "limon", "wheat", "lightblue", "lightpink", "deepolive", "red"]

        if domains:

            current_coordinates = (0, 0)
            for domain, coordinates in domains.items():

                c1 = set(range(current_coordinates[0], current_coordinates[1]))
                c2 = set(range(coordinates[0], coordinates[1]))

                # domain overlaps too much with the previous domain (33%) -> dont paint
                if len(c1 & c2) > len(c2) / 3:
                    continue

                # paint that domain
                else:
                    color = colors.pop(0)
                    f.write("color " + color + ", resi " + str(coordinates[0]) + "-" + str(coordinates[1]) + "\n")
                    middle = round(int((coordinates[1] - coordinates[0]) / 2 + coordinates[0]))
                    f.write('label (resi ' + str(middle) + ' and name CA), "%s" % ("' + str(domain) + '")\n')

                current_coordinates = coordinates

        # paint user residues if indicated
        if all_genes_dict:
            dup, label1, label2, label3 = all_genes_dict[gene]  # dup variable not used here
            # make sure label is present, end is bigger than start, end is within the protein length
            # do not allow start larger than end; do not allow end larger than total length

            if all(label1) \
                and (int(label1[1]) <= int(label1[2])) \
                and (int(label1[2]) <= int(list(matched_adaptive_sites_ref)[-1])):
                f.write("color " + "grey90" + ", resi " + label1[1] + "-" + label1[2] + "\n")
                middle = round((int(label1[2]) - int(label1[1])) / 2 + int(label1[1]))
                f.write('label (resi ' + str(middle) + ' and name CA), "%s" % ("' + label1[0] + '")\n')


            if all(label2) \
                and (int(label2[1]) <= int(label2[2])) \
                and (int(label2[2]) <= int(list(matched_adaptive_sites_ref)[-1])):
                f.write("color " + "yellow" + ", resi " + label2[1] + "-" + label2[2] + "\n")
                middle = round((int(label2[2]) - int(label2[1])) / 2 + int(label2[1]))
                f.write('label (resi ' + str(middle) + ' and name CA), "%s" % ("' + label2[0] + '")\n')

            if all(label3) \
                and (int(label3[1]) <= int(label3[2])) \
                and (int(label3[2]) <= int(list(matched_adaptive_sites_ref)[-1])):
                f.write("color " + "sand" + ", resi " + label3[1] + "-" + label3[2] + "\n")
                middle = round((int(label3[2]) - int(label3[1])) / 2 + int(label3[1]))
                f.write('label (resi ' + str(middle) + ' and name CA), "%s" % ("' + label3[0] + '")\n')

        # PyMOL command to color adaptive residues
        for site, features in matched_adaptive_sites_ref.items():

            # residues missing from PAML analysis
            if float(features[3]) == 0.00:
                residue = features[0] + str(site)
                f.write("select " + residue + ", resi " + str(site) + "\n")
                f.write("color gray40, " + residue + "\n")

            if float(features[2]) >= 0.90:
                residue = features[0] + str(site)
                if paint_sites:
                    f.write("select " + residue + ", resi " + str(site) + "\n")
                    f.write("color magenta, " + residue + "\n")
                    f.write("show sticks, " + residue + "\n")
                    f.write('label (resi ' + str(site) + ' and name CA), "%s" % ("' + residue + '")\n')

        # special case, first residue adaptive
        if float(matched_adaptive_sites_ref["1"][2]) >= 0.90:
            site = [site for site, features in matched_adaptive_sites_ref.items() if site == "1"][0]
            residue = matched_adaptive_sites_ref[site][0] + str(site)
            f.write('label (first (polymer and name CA)), "(%s; ' + residue + ')"%("N-term")\n')

        if float(matched_adaptive_sites_ref["1"][2]) < 0.90:
            f.write('label (first (polymer and name CA)), "(%s)"%("N-term")\n')

        # special case, last residue adaptive
        if float(matched_adaptive_sites_ref[str(len(matched_adaptive_sites_ref))][2]) >= 0.90:
            site = [site for site, features in matched_adaptive_sites_ref.items() if
                    site == str(len(matched_adaptive_sites_ref))][0]
            residue = matched_adaptive_sites_ref[site][0] + str(site)
            f.write('label (last (polymer and name CA)), "(%s; ' + residue + ')"%("C-term")\n')

        if float(matched_adaptive_sites_ref[str(len(matched_adaptive_sites_ref))][2]) < 0.90:
            f.write('label (last (polymer and name CA)), "(%s)"%("C-term")\n')

        # PyMOL command to set label size
        f.write("set label_size, 20\n")

        # PyMOL command to set label positions
        f.write("set label_position, (2,2,2)\n")

        # PyMOL command to set background color for saved output file
        f.write("set ray_opaque_background, on\n")

        # PyMOL command to save as figure (NO LICENSE PRINTS A NO LICENSE ON IMAGE)
        f.write("png " + protein_structure_path + gene + "_" + ref_species + ".png, width=12cm, height=8cm, "
                                                                             "dpi=300, ray=1\n")

        # PyMOL command to save the session
        f.write("save " + protein_structure_path + gene + "_" + ref_species + ".pse")

    return True

