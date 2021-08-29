#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 20:35:28 2021

@author: damian

Generates a PyMOL script and runs it internally. Saves the PyMOL session per protein.

"""

from ast import literal_eval
from Bio import SeqIO
import shutil
import subprocess
import os


# wdir = os.getcwd() + "/"
# protein = "Haus1"
# original_species = "Mm"
# result_path = wdir + "Results-06-13-2021-23-37/"

# LAST RESIDUE IS NOT MARKED IN PYMOL MODEL IF SCOREING (C-term label interferes?)
# Done but NOT TESTED YET


def check_structure(wdir, original_species, protein):
    """Checks presence of a structure prediction model for a given protein"""

    structure_model_path = wdir + "Structures/" + protein + "_" + original_species
    model_file_list = os.listdir(structure_model_path)
    unwanted_files = [file for file in model_file_list if file.startswith(".") or file.endswith(".png")]

    # check for unwanted files (hidden and png files) and remove them
    if unwanted_files:
        for file in unwanted_files:
            try:
                os.remove(structure_model_path + "/" + file)
            except FileNotFoundError:
                print("FileNotFoundError was triggered for: %s" % file)
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

    # there is exactly one pdb model
    #if len(model_file_list) == 1 and model_file_list[0].endswith(".pdb"):
    #    return True

    # there is more than one pdb model (not allowed)
    #if len(model_file_list) > 1:
    #    print("There is more than one structure prediction model for: %s -> skipping PyMOL for this protein" % protein)
    #    return False

    # there is no pdb model (not allowed)
    #if len(model_file_list) == 0:
    #    print("There is no structure prediction model for: %s -> skipping PyMOL for this protein" % protein)
    #    return False


def run_pymol(wdir, original_species, result_path, protein, proteins_under_positive_selection, offset):
    """Runs PyMOL with overlaid adaptive sites"""

    structures_path = wdir + "Structures"
    protein_path = structures_path + "/" + protein + "_" + original_species + "/"
    final_model_name = protein + "_" + original_species + ".pse"

    # get dict of adaptives sites matched to input sequence
    with open(result_path + "all_matched_adaptive_sites_original.txt", "r") as f:
        dictionary = literal_eval(f.read())

    # obtain a pymol script based on the model
    if not get_pymol_script(wdir, original_species, dictionary, protein, protein_path, proteins_under_positive_selection, offset):
        return False

    # run that script in pymol without triggering external GUI (-cq) -> DOES NOT WORK IN PYCHARM?
    pymol_command = "pymol -cq structure_overlay.pml"
    stderr, stdout = subprocess.Popen(pymol_command, shell=True, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).communicate()
    # move and overwrite if "structure_overlay.pml" exists in Structure folder for the protein
    shutil.move(os.path.join(wdir, "structure_overlay.pml"), os.path.join(protein_path, "structure_overlay.pml"))
    # move the model with overlayed residues into Results folder
    shutil.move(protein_path + final_model_name, result_path + final_model_name)

    return True


def get_pymol_script(wdir, original_species, dictionary, protein, protein_path, proteins_under_positive_selection, offset):
    """Gets a PyMOL script that will be passed into PyMOL automatically"""

    paint_sites = False

    # paint sites of proteins likely under positive selection
    if protein in proteins_under_positive_selection:
        paint_sites = True

    matched_adaptive_sites_original = dictionary[protein]
    structure_prediction_path = wdir + "Structures/" + protein + "_" + original_species

    if offset is None:
        offset = 0

    if len(os.listdir(structure_prediction_path)) == 0:
        print("\nNo structure predicion model is present for: %s -> cannot run PyMOL" % protein)
        return False
    #if len(os.listdir(structure_prediction_path)) > 1:
    #    print("\nMore than one structure predicion model is present for: %s -> cannot run PyMOL" % protein)
    #    return False

    with open("structure_overlay.pml", "w") as f:

        # PyMOL command to load model
        structure_model_path = wdir + "Structures/" + protein + "_" + original_species

        # select the model (check_structure function removed all hidden files)
        model = os.listdir(structure_model_path)[0]

        # start pymol script by loading the model
        f.write("load " + protein_path + model + "\n")

        # PyMOL command to color all resiues
        f.write("color cyan\n")

        # reindex all residues in the structure based on sequence used as a model structure
        f.write("alter (all), resi=str(int(resi)+" + str(offset) + ")\n")
        f.write("sort\n")
        f.write("rebuild\n")

        # PyMOL command to color adaptive residues
        for site, features in matched_adaptive_sites_original.items():

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

            if 0.90 > float(features[2]) >= 0.70:
                residue = features[0] + str(site)
                if paint_sites:
                    f.write("select " + residue + ", resi " + str(site) + "\n")
                    f.write("color lightblue, " + residue + "\n")
                    f.write("show sticks, " + residue + "\n")
                    f.write('label (resi '+ str(site) +' and name CA), "%s" % ("'+ residue +'")\n')

        # special case, first residue adaptive
        if float(matched_adaptive_sites_original["1"][2]) >= 0.70:
            site = [site for site, features in matched_adaptive_sites_original.items() if site == "1"][0]
            residue = matched_adaptive_sites_original[site][0] + str(site)
            f.write('label (first (polymer and name CA)), "(%s; ' + residue + ')"%("N-term")\n')

        if float(matched_adaptive_sites_original["1"][2]) < 0.70:
            f.write('label (first (polymer and name CA)), "(%s)"%("N-term")\n')

        # special case, last residue adaptive
        if float(matched_adaptive_sites_original[str(len(matched_adaptive_sites_original))][2]) >= 0.70:
            site = [site for site, features in matched_adaptive_sites_original.items() if
                    site == str(len(matched_adaptive_sites_original))][0]
            residue = matched_adaptive_sites_original[site][0] + str(site)
            f.write('label (last (polymer and name CA)), "(%s; ' + residue + ')"%("C-term")\n')

        if float(matched_adaptive_sites_original[str(len(matched_adaptive_sites_original))][2]) < 0.70:
            f.write('label (last (polymer and name CA)), "(%s)"%("C-term")\n')

        # PyMOL command to set label size
        f.write("set label_size, 20\n")

        # PyMOL command to set label positions
        f.write("set label_position, (2,2,2)\n")

        # PyMOL command to set background color for saved output file
        f.write("set ray_opaque_background, on\n")

        # PyMOL command to save as figure (NO LICENSE PRINTS A NO LICENSE ON IMAGE)
        f.write("png " + protein_path + protein + "_" + original_species + ".png, width=12cm, height=8cm, dpi=300, ray=1\n")

        # PyMOL command to save the session
        f.write("save " + protein_path + protein + "_" + original_species + ".pse")

    return True


"""

# THIS IS NOT NEEDED ANYMORE:
def compare_model_with_input(wdir, original_species, protein, model_seq):
    
    # silence warnings from Biopython about missing header in model (has to follow import)
    #import warnings
    #from Bio import BiopythonWarning
    #warnings.simplefilter('ignore', BiopythonWarning)
    
    # get input protein sequence used for blast
    blast_input_path = wdir + "Blast_input/"
    input_seq_filename = protein + "_" + original_species + "_protein.fasta"
    with open(blast_input_path + input_seq_filename) as f:
        input_seq = f.readlines()[1]

    # get protein sequence of the model
    #structure_model_path = wdir + "Structures/" + protein + "_" + original_species
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
def check_all_structures(wdir, original_species):
    Checks presence of structure prediction models for all proteins
   
    all_proteins = [protein.rstrip("\n") for protein in open(wdir + "proteins.txt", "r").readlines()]
    missing_structures = [check_structure(wdir, original_species, protein) for protein in all_proteins]
    missing_structures_final = [structure for structure in missing_structures if structure is not None]
    
    if missing_structures_final:
        print("...WARNING... (FREEDA) I did not find clear structure prediction models for: %s" % missing_structures_final)
        print("...WARNING... (FREEDA) I cannot overlay adaptive sites for these proteins\n")

    return missing_structures_final

"""
