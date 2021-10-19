#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 20:35:28 2021

@author: damian

Generates a PyMOL script and runs it internally. Saves the PyMOL session per protein.

"""

from freeda import fasta_reader
from ast import literal_eval
import shutil
import subprocess
import os
import requests
import simplejson.errors


# wdir = os.getcwd() + "/"
# protein = "Haus1"
# ref_species = "Mm"
# result_path = wdir + "Results-06-13-2021-23-37/"

# LAST RESIDUE IS NOT MARKED IN PYMOL MODEL IF SCORING (C-term label interferes?)
# Done but NOT TESTED YET



# 2021-09-15
# Brian Akins

# METHOD: get_interpro(uniprot_id), interacts with the InterPro REST API.
# INPUT:  A string containing the Uniprot ID of a target protein.
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
    #print(response.json())
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


def get_domains(wdir, ref_species, protein):
    """Based on uniprot id gets domain layout from Interpro API : name and coordinates"""

    structure_filepath = wdir + "Structures/" + protein + "_" + ref_species

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

        # no Interpro data for this protein
        except simplejson.errors.JSONDecodeError:
            return domains
        # requested url wasnt pulled correctly
        except KeyError:
            count += 1

    return domains


def check_structure(wdir, ref_species, protein):
    """Checks presence of a structure prediction model for a given protein"""

    structure_model_path = wdir + "Structures/" + protein + "_" + ref_species
    model_file_list = os.listdir(structure_model_path)
    unwanted_files = [file for file in model_file_list if file.startswith(".") or file.endswith(".png")]

    # check for unwanted files (hidden and png files) and remove them
    if unwanted_files:
        for file in unwanted_files:
            try:
                os.remove(structure_model_path + "/" + file)
            except FileNotFoundError:
                #print("FileNotFoundError was triggered for: %s" % file)
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


def run_pymol(wdir, ref_species, result_path, protein, proteins_under_positive_selection, offset):
    """Runs PyMOL with overlaid adaptive sites"""

    structures_path = wdir + "Structures"
    protein_path = structures_path + "/" + protein + "_" + ref_species + "/"
    final_model_name = protein + "_" + ref_species + ".pse"

    # get dict of adaptives sites matched to input sequence
    with open(result_path + "all_matched_adaptive_sites_ref.txt", "r") as f:
        dictionary = literal_eval(f.read())

    # get domain names and coordinates from Interpro API
    domains = get_domains(wdir, ref_species, protein)

    # obtain a pymol script based on the model
    if not get_pymol_script(wdir, ref_species, dictionary, protein,
                            protein_path, proteins_under_positive_selection, offset, domains):
        return False

    # run that script in pymol without triggering external GUI (-cq) -> DOES NOT WORK IN PYCHARM?
    pymol_command = "pymol -cq structure_overlay.pml"
    stderr, stdout = subprocess.Popen(pymol_command, shell=True, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).communicate()
    # move and overwrite if "structure_overlay.pml" exists in Structure folder for the protein
    shutil.move(os.path.join(wdir, "structure_overlay.pml"), os.path.join(protein_path, "structure_overlay.pml"))
    # move the model with overlaid residues into Results folder
    shutil.move(protein_path + final_model_name, result_path + final_model_name)

    return True


def get_pymol_script(wdir, ref_species, dictionary, protein,
                     protein_path, proteins_under_positive_selection, offset, domains):
    """Gets a PyMOL script that will be passed into PyMOL automatically"""

    paint_sites = False

    # paint sites of proteins likely under positive selection
    if protein in proteins_under_positive_selection:
        paint_sites = True

    matched_adaptive_sites_ref = dictionary[protein]
    structure_prediction_path = wdir + "Structures/" + protein + "_" + ref_species

    if offset is None:
        offset = 0

    if len(os.listdir(structure_prediction_path)) == 0:
        print("\nNo structure predicion model is present for: %s -> cannot run PyMOL" % protein)
        return False

    with open("structure_overlay.pml", "w") as f:

        # PyMOL command to load model
        structure_model_path = wdir + "Structures/" + protein + "_" + ref_species

        # select the model (check_structure function removed all hidden files)
        model = [file for file in os.listdir(structure_model_path) if file.endswith(".pdb")][0]

        # start pymol script by loading the model
        f.write("load " + protein_path + model + "\n")

        # PyMOL command to color all residues
        f.write("color cyan\n")

        # reindex all residues in the structure based on sequence used as a model structure
        f.write("alter (all), resi=str(int(resi)+" + str(offset) + ")\n")
        f.write("sort\n")
        f.write("rebuild\n")

        # color domains
        colors = ["yellow", "orange", "marine", "limon", "wheat", "lightblue", "lightpink", "deepolive", "red"]

        if domains:

            current_coordinates = (0, 0)
            for domain, coordinates in domains.items():

                c1 = set(range(current_coordinates[0], current_coordinates[1]))
                c2 = set(range(coordinates[0], coordinates[1]))

                # domain overlaps too much with the previous domain -> dont paint
                if len(c1 & c2) > len(c2) / 2:
                    continue

                # paint that domain
                else:
                    color = colors.pop(0)
                    f.write("color " + color + ", resi " + str(coordinates[0]) + "-" + str(coordinates[1]) + "\n")
                    middle = int((coordinates[1] - coordinates[0]) / 2 + coordinates[0])
                    f.write('label (resi ' + str(middle) + ' and name CA), "%s" % ("' + str(domain) + '")\n')

                current_coordinates = coordinates

        # PyMOL command to color adaptive residues
        for site, features in matched_adaptive_sites_ref.items():

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

            #if 0.90 > float(features[2]) >= 0.70:
            #    residue = features[0] + str(site)
            #    if paint_sites:
            #        f.write("select " + residue + ", resi " + str(site) + "\n")
            #        f.write("color lightblue, " + residue + "\n")
            #        f.write("show sticks, " + residue + "\n")
            #        f.write('label (resi '+ str(site) +' and name CA), "%s" % ("'+ residue +'")\n')

        # special case, first residue adaptive
        if float(matched_adaptive_sites_ref["1"][2]) >= 0.70:
            site = [site for site, features in matched_adaptive_sites_ref.items() if site == "1"][0]
            residue = matched_adaptive_sites_ref[site][0] + str(site)
            f.write('label (first (polymer and name CA)), "(%s; ' + residue + ')"%("N-term")\n')

        if float(matched_adaptive_sites_ref["1"][2]) < 0.70:
            f.write('label (first (polymer and name CA)), "(%s)"%("N-term")\n')

        # special case, last residue adaptive
        if float(matched_adaptive_sites_ref[str(len(matched_adaptive_sites_ref))][2]) >= 0.70:
            site = [site for site, features in matched_adaptive_sites_ref.items() if
                    site == str(len(matched_adaptive_sites_ref))][0]
            residue = matched_adaptive_sites_ref[site][0] + str(site)
            f.write('label (last (polymer and name CA)), "(%s; ' + residue + ')"%("C-term")\n')

        if float(matched_adaptive_sites_ref[str(len(matched_adaptive_sites_ref))][2]) < 0.70:
            f.write('label (last (polymer and name CA)), "(%s)"%("C-term")\n')

        # PyMOL command to set label size
        f.write("set label_size, 20\n")

        # PyMOL command to set label positions
        f.write("set label_position, (2,2,2)\n")

        # PyMOL command to set background color for saved output file
        f.write("set ray_opaque_background, on\n")

        # PyMOL command to save as figure (NO LICENSE PRINTS A NO LICENSE ON IMAGE)
        f.write("png " + protein_path + protein + "_" + ref_species + ".png, width=12cm, height=8cm, dpi=300, ray=1\n")

        # PyMOL command to save the session
        f.write("save " + protein_path + protein + "_" + ref_species + ".pse")

    return True





#cenpt_uniprot = 'Q3TJM4'
#cenpo_uniprot = 'Q9BU64'
#test_cenpt = protein_domains(cenpt_uniprot)
#test_cenpo = protein_domains(cenpo_uniprot)
#print(test_cenpt)
#print(test_cenpo)

"""

# THIS IS NOT NEEDED ANYMORE:
def compare_model_with_input(wdir, ref_species, protein, model_seq):
    
    # silence warnings from Biopython about missing header in model (has to follow import)
    #import warnings
    #from Bio import BiopythonWarning
    #warnings.simplefilter('ignore', BiopythonWarning)
    
    # get input protein sequence used for blast
    blast_input_path = wdir + "Blast_input/"
    input_seq_filename = protein + "_" + ref_species + "_protein.fasta"
    with open(blast_input_path + input_seq_filename) as f:
        input_seq = f.readlines()[1]

    # get protein sequence of the model
    #structure_model_path = wdir + "Structures/" + protein + "_" + ref_species
    #model_filename = [model for model in os.listdir(structure_model_path) if model.endswith(".pdb")][0]
    #model_path = structure_model_path + "/" + model_filename
    
    #with open(model_path, 'r') as pdb_file:
    #    for record in SeqIO.parse(pdb_file, 'pdb-atom'):
    #        model_seq = record.seq

    # compare sequences and do not allow model overlay if sequences differ
    if input_seq != model_seq:
        print("\n...WARNING... Protein sequence generated DOES NOT match the model for: %s\n" \
                                              "-> cannot run PyMOL\n" % protein)
        print("Input sequence:\n%s\n" % input_seq)
        print("Model sequence:\n%s\n" % model_seq)
        return
    
    else:
        print("\nProtein sequence generated matches the model for protein: %s" % protein)
        return


# THIS IS PROBABLY NOT NEEDED:
def check_all_structures(wdir, ref_species):
    Checks presence of structure prediction models for all proteins
   
    all_proteins = [protein.rstrip("\n") for protein in open(wdir + "proteins.txt", "r").readlines()]
    missing_structures = [check_structure(wdir, ref_species, protein) for protein in all_proteins]
    missing_structures_final = [structure for structure in missing_structures if structure is not None]
    
    if missing_structures_final:
        print("...WARNING... (FREEDA) I did not find clear structure prediction models for: %s" % missing_structures_final)
        print("...WARNING... (FREEDA) I cannot overlay adaptive sites for these proteins\n")

    return missing_structures_final

"""
