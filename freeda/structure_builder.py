#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 20:35:28 2021

@author: damian

Generates a PyMOL script and runs it internally. Saves the PyMOL session per protein.

"""

from ast import literal_eval
import shutil
import subprocess
import os

#wdir = os.getcwd() + "/"
#result_path = wdir + "Results-06-13-2021-23-37/"

"""
from Bio import PDB
from Bio.PDB.MMCIFParser import MMCIFParser

pdbl = PDB.PDBList()
pdbl.retrieve_pdb_file("4o9f")

parser = MMCIFParser()
structure = parser.get_structure("4o9f", "4o9f.cif")




structure_id = "4O9F".lower()
filename = "pdb4o9f.ent"
structure = parser.get_structure(structure_id, filename)

"""


def run_pymol(wdir, original_species, result_path, offset):
    
    structures_path = wdir + "Structures"
    
    with open(result_path + "all_matched_adaptive_sites_original.txt", "r") as f:
        dictionary = literal_eval(f.read())
    
    with open("proteins.txt", "r") as f:
       file = f.readlines()
        
       for line in file:
           protein = line.rstrip("\n")
           protein_path = structures_path + "/" + protein + "_" + original_species + "/"
        
           # obtain a pymol script based on the model
           # protein_path in Structures folder MUST EXIST !!!
           get_pymol_script(wdir, result_path, dictionary, protein, protein_path, offset)
       
           # run that script in pymol without triggering external GUI (-cq)
           pymol_command = "pymol -cq structure_overlay.pml"
           stderr, stdout = subprocess.Popen(pymol_command, shell=True, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
           # move and overwrite if "structure_overlay.pml" exists in Structure folder for the protein
           shutil.move(os.path.join(wdir, "structure_overlay.pml"), os.path.join(protein_path, "structure_overlay.pml"))
        
        
def get_pymol_script(wdir, result_path, dictionary, protein, protein_path, offset):
    
    matched_adaptive_sites_original = dictionary[protein]
    
    with open("structure_overlay.pml", "w") as f:
        
        # PyMOL command to load model
        f.write("load " + protein_path + "model1.pdb\n")
    
        # PyMOL command to color all resiues
        f.write("color cyan\n")
        
        # reindex all residues in the structure based on sequence used as a model structure
        f.write("alter (all), resi=str(int(resi)+" + str(offset) + ")\n")
        f.write("sort\n")
        f.write("rebuild\n")
        
        # PyMOL command to color adaptive residues
        for site, features in matched_adaptive_sites_original.items():
            if float(features[2]) >= 0.90:
                residue = features[0] + str(site)
                f.write("select " + residue + ", resi " + str(site) + "\n")
                f.write("color magenta, " + residue + "\n")
                f.write("show sticks, " + residue + "\n")
                f.write('label (resi '+ str(site) +' and name CA), "%s" % ("'+ residue +'")\n')
            
            if 0.90 > float(features[2]) >= 0.75:
                residue = features[0] + str(site)
                f.write("select " + residue + ", resi " + str(site) + "\n")
                f.write("color white, " + residue + "\n")
                f.write("show sticks, " + residue + "\n")
                f.write('label (resi '+ str(site) +' and name CA), "%s" % ("'+ residue +'")\n')
        
        # PyMOL comand to mark N-term and C-term
        f.write('label (first (polymer and name CA)), "(%s)"%("N-term")\n')
        f.write('label (last (polymer and name CA)), "(%s)"%("C-term")\n')
            
        # PyMOL command to set label size
        f.write("set label_size, 20\n")
            
        # PyMOL command to set label positions
        f.write("set label_position, (2,2,2)\n")
            
        # PyMOL command to set background color for saved output file
        f.write("set ray_opaque_background, on\n")
        
        # PyMOL command to save as figure (NO LICENSE PRINTS A NO LICENSE ON IMAGE)
        f.write("png " + protein_path + protein + ".png, width=12cm, height=8cm, dpi=300, ray=1\n")
            
        # PyMOL command to save the session
        f.write("save " + protein_path + protein + ".pse")
    