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

#wdir = os.getcwd() + "/"
#result_path = wdir + "Results-06-13-2021-23-37/"

def run_pymol(wdir, result_path):
    
    structures_path = wdir + "Structures"
    
    with open(result_path + "all_matched_adaptive_sites_original.txt", "r") as f:
        dictionary = literal_eval(f.read())
    
    with open("proteins.txt", "r") as f:
       file = f.readlines()
        
       for line in file:
           protein = line.rstrip("\n")
           protein_path = structures_path + "/" + protein + "/"
        
           # obtain a pymol script based on the model
           # protein_path in Structures folder MUST EXIST !!!
           get_pymol_script(wdir, result_path, dictionary, protein, protein_path)
       
           # run that script in pymol without triggering external GUI (-cq)
           pymol_command = "pymol -cq structure_overlay.pml"
           stderr, stdout = subprocess.Popen(pymol_command, shell=True, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
           shutil.move("structure_overlay.pml", protein_path)
        
def get_pymol_script(wdir, result_path, dictionary, protein, protein_path):
    
    # STILL WORKING ON THAT FUNCITON -> outputs a PyMOL script into Data folder
    # to copy and paste into PyMOL
    
    matched_adaptive_sites_original = dictionary[protein]
    
    with open("structure_overlay.pml", "w") as f:
        
        # PyMOL command to load model
        f.write("load " + protein_path + "model1.pdb\n")
    
        # PyMOL command to color all resiues
        f.write("color cyan\n")
            
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
    