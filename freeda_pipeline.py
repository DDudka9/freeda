#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:36:16 2021

@author: damian

Before testing make sure that code can communicate witg files:
    1) Make sure that genomes are called from a separate folder
    2) Check generate folders module
    3) Change directories of Output files for blast

"""
# 06/17/2021
# IMPORTANT -> check matches generator module -> does it sort matches before chopping into smaller?
# check on nomleu3 genome CD46 chr5_rev -> compare to gorGor6 chr1_6_for
# try changing min threshold for matches from 40 to 60
# I think the problem is in sorting matches and the fact that two domains in CD46 are quite repeated
# so broken large contigs do not align well and distant exons 10-13 are not found on later contigs
# this is speific to exons remote by over 13kb
# consider increasing suffix and prefix to 15kb? But this will not help with repeated domains
# even if exon 10 will be found then exon 11 will not be aligned properly 
# WRITE WARNING FOR PROTEINS WITH REPEATED DOMAINS (ex. TACC3, CD46)

# 03_28_2021
# Cenpt Caroli chr19__rev somehow sees only 12 exons (not 13) and final
# exon count is set as 12 not 13 (even if the right contig picked chr8__rev has 13 exons)
# its clear that cloner depends on the last contig analysed which is wrong
# but mainly its a failure of the find_exons function
# try to get exons dictionary once and then use try/except to avoid KeyError
# when cloning exons (and potentially having different number)
# SOLUTION: find_exons doesnt output final exon count anymore (FIXED)
# and find_exons function requires more fixes (in progress)

# Cenpt runs well with the packaged FREEDA (almost the same M2a and M8)

# Generating names for folders from genomes names is defective
# I think its cose "lstrip" that strips "Acomys_Rus -> Acomys_Ru" and "Neomys_lepida" -> "Neomys_lepid"

# VERSION CONTROL:
    # 1) Go to /Volumes/DamianEx/FREEDA_v01/src (or ...src/FREEDA)
    # 2) In IPython console type:
    # 3) !git add "main.py" (or whatever the module you changed)
    # 4) !git commit -m "My commit"
    # 5) !git push origin main (need to make another brunch except main)

# !git commit -m 'Initial commit'
# !git add forgotten_file
# !git commit --amend

# !git branch <newbranchname> (create a new branch)
# !git checkout -b <newbranchname> (create and got to new branch)
# !git checkout main (go to a branch)
# !git checkout - (go back to pewviously checkedout branch)
# always have clear staging area (all commited) before switching branches
# !git merge <newbranchname> (merge the new branch with main)
# !git branch -d <newbranchname> (delete the new branch cose main points already at the same place)
# dont modify the same place in two different branches!!!!

# Thre is a problem in recognising the first exon in Terf2 Apodemus sylv contig LIPJ01013670.1__rev
# leading hypothesis is that the exon does not fall into any of the "callings"
# dissect it and run separately with printing functions


# Define a modeule for tweaking parameters ex.
#   - duplication restriction (switches on the duplication score)
#   - blast threshold
#   - homology threshold
#   - synteny threshold
#   - coverage threshold
#   - non_ACGT corrector (to mirror CDS position)

# I ADDED DELE&ION TO INSERTION CHECK
# I ALSO DISABLED FRAMESHIFT CHECK in cds_cloner-py 
# UNLESS ITS A DUPLICATED EXON (primates seem to have legit indels)

# HAUS2 ponAbe3 -> chr2A has many more exons cose cds != genpmic in alignment (pieces)
# enforce that exon needs to not onlx != "-" but also be equal to ake exon
# also HAUS3 genomic sequence from ensembl is missing the last exon !?!? (download again)
# FIXED -> made a mistake providing genomic sequence (shorter)


# Something weird about Bub1 -> lots of >0.90 sites but M7 higher than M8

# Cloner module needs revision to get hamming distance duplication comparison compare
# the actual duplicated exons and not only the number of exon they carry
# test on Aurkc Ap

# Last bp in Anapc16 is plotted as red on graph -> why?
# Gblocks seems to not like the last codon for some reason (last after STOP removal)

# Showcase STOP_remover function using Cenol -> Pd has an early STOP
# which is followed by a mutation in the original STOP suggesting a rescue that way
# in the same time Mc also has a mutation in the STOP but doesnt have an early STOP rescue
# suggesting that there is another later STOP that FREEDA does not capture (WHICH THERE ISNT!!!)

# Single non_ACGT bases currently lead to whol exon loss
# SOLUTION: THINK ABOUT FLIPPING non_ACGT INTO CORRESPONDING CDS POSITION (conservative)
# this could save these exons!

# THERE IS AN ISSUE WITH: if earlier STOP present in other species then
# original species gets translated normally and final_original_dict is +1
# which leads to ValueError in get_omegas function 
# SOLUTION: enabled the STOP_remover function( FIXED?)
# TO FIX: dashes in the MAFFT alignment (need to remove these positions before
# counting codons -> test on Haus8)

# Think about giving 0.5 RETRO score value for Retro at one end and introny at another
# to allow hamming distance function to kick in (not sure about this)


# THERE IS AN ERROR IN HAUS8 CORRECTION function -> not same lengths?
# check the print screen


# How come there is no contig with Cenpt coming from Spicilegus genome?
# Make FREEDA log file more readable (indentations)


from freeda import input_extractor
from freeda import tblastn
from freeda import exon_extractor
from freeda import paml_launcher
from freeda import paml_visualizer
from freeda import structure_builder

import os

def freeda_pipeline(original_species=None, t=None, mafft_path=None):
    
    # this needs to be a path to the "Data" folder
    wdir = os.getcwd() + "/"
    
    # reference species sequences: protein seq, cds, exons, gene (ex. Mus musculus)
    if original_species == None:
        original_species = "Mm"
    
    # initial percent identity threshold for blast matches analysis
    if t == None:
        t = 30
    
    # this needs to be a path where MAFFT aligner is installed
    if mafft_path == None:
        mafft_path = "/Users/damian/anaconda3/bin/mafft"
    
    user_input1 = None
    user_input2 = None
    user_input3 = None
    user_input4 = None
    structure_overlay = True
    
    structure_model_present = check_structure(wdir, original_species)
    
    
######## GET USER INPUT ########
    
    
    while user_input1 != "y" and user_input1 != "n":
        user_input1 = input("(FREEDA) Should I run blast? (y / n)\n").lower()
        if user_input1.lower() != "y" and user_input1.lower() != "n":
            print("Please answer y or n\n")
        
    while user_input2 != "y" and user_input2 != "n":
        user_input2 = input("(FREEDA) Should I find exons? (y / n)\n").lower()
        if user_input2.lower() != "y" and user_input2.lower() != "n":
            print("Please answer y or n\n")
        
    while user_input3 != "y" and user_input3 != "n":
        user_input3 = input("(FREEDA) Should I run PAML? (y / n)\n").lower()
        if user_input3.lower() != "y" and user_input3.lower() != "n":
            print("Please answer y or n\n")
    
    while user_input4 != "y" and user_input4 != "n":
        user_input4 = input("(FREEDA) Should I run PyMOL? (y / n)\n").lower()
        if user_input4.lower() != "y" and user_input4.lower() != "n":
            print("Please answer y or n\n")
    
    # make sure all structure prediction models are present
    if user_input4 == "y" and not structure_model_present:
        structure_overlay = False
        print("...WARNING... (FREEDA) I cannot overlay adaptive sites until all proteins have structure predictions\n")
        print("...WARNING... (FREEDA) I will run the pipeline skipping PyMOL\n")
        print("...WARNING... (FREEDA) You can run PyMOL later on based on these results\n")
    
    if user_input4 == "y" and structure_model_present:
        offset = input("(FREEDA) What is the offset used for modelling? (number)\n")
    
    
######## GET ALL INPUT DATA  ########

    if user_input1 == "y":

        # generate a reference Genome object
        reference_genome_present, ensembl, original_species, reference_genomes_path, reference_genome_name, \
            reference_genome_contigs_dict, biotype = input_extractor.generate_reference_genome_object(wdir, original_species)
            
        if reference_genome_present == True:
            all_proteins = open(wdir + "proteins.txt", "r").readlines()
            for protein in all_proteins:
                input_extractor.extract_input(wdir, original_species, reference_genome_name, reference_genomes_path, 
                        reference_genome_contigs_dict, ensembl, biotype, protein.rstrip("\n"))
            print("\nInput data have been generated -> checking genome blast databases...\n")
        
        if reference_genome_present == False:
            print("\n(FREEDA) I couldnt find the reference genome" \
              "\n   Make sure you downloaded it into .../Data/Reference_genomes from https://www.ncbi.nlm.nih.gov/assembly (mouse: GCA_000001635.8; human: GCA_000001405.28) "\
                  "-> exiting the pipeline now...")
            return
    
    
######## RUN BLAST ########
    

    if user_input1 == "y":
        blast_path = tblastn.run_blast(wdir, original_species)
        if blast_path == None:
            print("\nBlast database build failed for at least one genome" \
                  "\n   Make sure you downloaded all genomes -> exiting the pipeline now...")
            return
        
    else:
        blast_path = wdir + "Blast_output/"



    # ADD A CHECK FOR BLAST OUTPUT IN CASE IT WAS TEMPERED WITH OR SPECIES WERE REMOVED



######## RUN EXON FINDING ########


    
    if user_input2 == "y":
        result_path = exon_extractor.analyse_blast_results(wdir, blast_path, \
                                        original_species, int(t), mafft_path)



######## RUN PAML ########



    if user_input3 == "y" and user_input2 == "n":
        nr_of_tries = 1
        while nr_of_tries <= 3:
            user_input4 = input("----->  Indicate folder with results (no slashes or quotes): ")
            result_path = wdir + user_input4 + "/"
            if os.path.isdir(result_path) == False:
                print("\n(FREEDA) I couldnt find your results folder (%s/3)" % nr_of_tries)
                nr_of_tries += 1
            else:
                nr_of_tries = float("inf")
                result_path = wdir + user_input4 + "/"
                proteins, nr_of_species_total_dict, PAML_logfile_name, day = paml_launcher.analyse_final_cds(wdir, 
                                                                original_species, result_path, mafft_path)
                paml_visualizer.analyse_PAML_results(wdir, result_path, 
                            proteins, nr_of_species_total_dict, original_species, PAML_logfile_name, day)
                
    if user_input3 == "y" and user_input2 == "y":
        proteins, nr_of_species_total_dict, PAML_logfile_name, day = paml_launcher.analyse_final_cds(wdir, 
                                                            original_species, result_path, mafft_path)
        paml_visualizer.analyse_PAML_results(wdir, result_path, 
                            proteins, nr_of_species_total_dict, original_species, PAML_logfile_name, day)
    
    
######## RUN PyMOL ########
    
    
    if user_input1 == "n" and user_input2 == "n" and user_input3 == "n" and user_input4 == "y" and structure_model_present:
        nr_of_tries = 1
        while nr_of_tries <= 3:
            user_input4 = input("----->  Indicate folder with results (no slashes or quotes): ")
            result_path = wdir + user_input4 + "/"
            if os.path.isdir(result_path) == False:
                print("\n(FREEDA) I couldnt find your results folder (%s/3)" % nr_of_tries)
                nr_of_tries += 1
            else:
                nr_of_tries = float("inf")
                result_path = wdir + user_input4 + "/"
                structure_builder.run_pymol(wdir, original_species, result_path, offset)
    
                # indicate that structure overlay was already run
                structure_overlay = False
                
    if structure_model_present and structure_overlay:
        structure_builder.run_pymol(wdir, original_species, result_path, offset)
        
    
    print("\nYou reached the end of FREEDA pipeline.")




def check_structure(wdir, original_species):
    
    with open("proteins.txt", "r") as f:
        file = f.readlines()
        
        for line in file:
            protein = line.rstrip("\n")
            structures_path = wdir + "Structures"
            model_path = structures_path + "/" + protein + "_" + original_species + "/model1.pdb"
        
            if os.path.isfile(model_path):
                structure_model_present = True
            else:
                structure_model_present = False
                print("\n...WARNING... Structure prediction for protein: %s DOES NOT EXIST" % protein)
    
    # do not return any script if even one model is not present
    return structure_model_present


if __name__ == '__freeda_pipeline__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Run FREEDA pipeline")
    parser.add_argument("-os", "--original_species", 
        help="specify reference organism (ex. Mm for Mus musculus)", type=str,
        required=False)
    parser.add_argument("-t", "--blast_threshold", 
        help="specify percentage identity threshold for blast (ex. 30; default)", type=int,
        required=True)
    parser.add_argument("-mafft", "--mafft_path", help="specify path to mafft aligner (ex. /Users/user/anaconda/bin/mafft)", type=str,
        required=True)
    
    args = parser.parse_args()
    freeda_pipeline(original_species=args.original_species, t=args.blast_threshold, mafft_path=args.mafft_path)

    
    
    
    
    
    
    
    
    