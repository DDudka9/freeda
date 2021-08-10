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

# How come "Contig too short to check C-term synteny 0bp aligned" for contig 81143 in genome11 Ap2m1
# SOLUTION: Probably connected to exon4 being a 6bp microexon and NOT deleted from exon input but
# There are 13 exons expected instead of 11 -> exon 12 is skipped for some reason; alignment looks good
# Also alignment of single exons from exon 7 is messed up (linux default file order problem again?)
# early STOP remover function worked well -> post trimming it was easier to align hence difference in "no_STOP" alignment length
# but since last 4 single exons were aligned poorly, the stop codons were missing/were displaced in other species

# GET RID OF MAFFT PATH -> not needed if mafft binary in the path (with other dependecies)
# DONE

# Try downloading all models for a species into freeda code file? Zipped? Getting models without unzipping?

# First check model sequence and length -> pick pyensembl transcripr based on that or longest (no user prompt)


from freeda import input_extractor
from freeda import tblastn
from freeda import exon_extractor
from freeda import paml_launcher
from freeda import paml_visualizer
from freeda import structure_builder
from json import dump
from ast import literal_eval
import os
import shutil


def freeda_pipeline(original_species=None, t=None):
    # current directory must be the "Data" folder

    wdir = os.getcwd() + "/"

    # reference species sequences: protein seq, cds, exons, gene (ex. Mus musculus)
    if original_species is None:
        original_species = "Mm"

    # initial percent identity threshold for blast matches analysis
    if t is None:
        t = 30

    user_input0 = None
    user_input1 = None
    user_input2 = None
    user_input3 = None
    user_input4 = None
    structure_prediction_matching = False
    structure_overlay_present = False

    # ----------------------------------------#
    ######## GET USER INPUT ########
    # ----------------------------------------#

    while user_input0 != "y" and user_input0 != "n":
        user_input0 = input("(FREEDA) Should I get input? (y / n)\n").lower()
        if user_input0.lower() != "y" and user_input0.lower() != "n":
            print("Please answer y or n\n")

    while user_input1 != "y" and user_input1 != "n":
        user_input1 = input("(FREEDA) Should I run blast? (y / n)\n").lower()
        if user_input1.lower() != "y" and user_input1.lower() != "n":
            print("Please answer y or n\n")

    while user_input2 != "y" and user_input2 != "n":
        user_input2 = input("(FREEDA) Should I find exons? (y / n)\n").lower()
        if user_input2.lower() != "y" and user_input2.lower() != "n":
            print("Please answer y or n\n")

    while user_input3 != "y" and user_input3 != "n":
        user_input3 = input("(FREEDA) Should I perform molecular evolution analysis (PAML)? (y / n)\n").lower()
        if user_input3.lower() != "y" and user_input3.lower() != "n":
            print("Please answer y or n\n")

    while user_input4 != "y" and user_input4 != "n":
        user_input4 = input(
            "(FREEDA) Should I overlay putative adaptive sites onto 3D structure (PyMOL)? (y / n)\n").lower()
        if user_input4.lower() != "y" and user_input4.lower() != "n":
            print("Please answer y or n\n")

    # if user_input4 == "y":
    #    offset = input("(FREEDA) What is the offset used for modelling? (e.g. type 1 if full length)\n")

    # ----------------------------------------#
    ######## GET ALL INPUT DATA  ########
    # ----------------------------------------#

    # record inputed data parameters for using in later stages of the pipeline
    input_dictionary = {}

    # User wants to generate input data
    if user_input0 == "y":

        # generate a reference Genome object
        reference_genome_name = input("(FREEDA) What is the name of the reference genome? (e.g. MUSCULUS_genome)\n")
        reference_genome_present, ensembl, original_species, reference_genomes_path, reference_genome_contigs_dict, \
        biotype = input_extractor.generate_reference_genome_object(wdir, original_species, str(reference_genome_name))

        # stop pipeline if the reference genome is absent
        if not reference_genome_present:
            print("\n...FATAL ERROR...: There is no reference genome detected -> exiting the pipeline now...\n"
                  "\n   Make sure you downloaded it into ../Data/Reference_genomes from "
                  " https://www.ncbi.nlm.nih.gov/assembly -> (mouse: GCA_000001635.8; human: GCA_000001405.28) -> "
                  "GenBank -> Genomic FASTA(.fna)")
            return

        all_proteins = [protein.rstrip("\n") for protein in open(wdir + "proteins.txt", "r").readlines()]
        for protein in all_proteins:
            if protein == "\n":
                continue
            print("\n----------- * %s * -----------" % protein)

            # get structure prediction model from AlphaFold
            input_dictionary[protein] = []
            possible_uniprot_ids = input_extractor.get_uniprot_id(wdir, original_species, protein)
            model_seq = input_extractor.fetch_structure_prediction(wdir, original_species, protein,
                                                                   possible_uniprot_ids)
            # input_dictionary[protein].append(str(model_seq))

            # get sequence input from ensembl
            input_correct, model_matches_input, microexon_present = input_extractor.extract_input(wdir,
                                                                                                  original_species,
                                                                                                  reference_genome_name,
                                                                                                  reference_genomes_path,
                                                                                                  reference_genome_contigs_dict,
                                                                                                  ensembl, biotype,
                                                                                                  protein, model_seq)

            if input_correct:
                print("\nInput data have been generated for protein: %s\n\n" % protein)

            if not input_correct:
                print(
                    "\n...FATAL ERROR...: Input data generation FAILED for protein: %s -> exiting the pipeline now...\n" % protein)
                return

            if model_matches_input and not microexon_present:
                # ast module requires a string
                input_dictionary[protein] = str(model_matches_input)

            if not model_matches_input:
                print(
                    "...WARNING...: Structure prediction for protein: %s DOES NOT have a match in available ensembl "
                    "database -> cannot run PyMOL\n" % protein)
                print("...WARNING...: Protein will be analyzed using PAML without 3D structure overlay\n")
                # ast module requires a string
                input_dictionary[protein] = str(model_matches_input)

            if microexon_present:
                model_matches_input = False
                print("...WARNING...: Sequence for: %s found in Ensembl contains a microexon\n" % protein)
                print("...WARNING...: Microexons are difficult to align and are removed -> cannot run PyMOL\n")
                # ast module requires a string
                input_dictionary[protein] = str(model_matches_input)

        # save the input dict into a file
        with open("most_recent_input_dictionary.txt", "w") as file:
            dump(input_dictionary, file)
        shutil.move(wdir + "most_recent_input_dictionary.txt", wdir + "Blast_input/most_recent_input_dictionary.txt")

        print("\nAll input data have been generated\n")

    # User doesnt want to generate input data
    if user_input0 == "n":

        # check if input data have been generated beforehand
        path_to_input_dict = wdir + "Blast_input/"
        if os.path.exists(path_to_input_dict + "most_recent_input_dictionary.txt"):
            with open(path_to_input_dict + "most_recent_input_dictionary.txt", "r") as f:
                # ast module requires a string
                input_dictionary = literal_eval(f.read())
        else:
            print("\n...FATAL ERROR...: You need to generate input data first -> exiting the pipeline now...")
            return

        # DO I NEED THIS? Models are always present but should not be used sometimes -> input_dictionary.txt:

        # check if all models are present
        # all_proteins = [protein.rstrip("\n") for protein in open(wdir + "proteins.txt", "r").readlines()]
        # missing_structures = [structure_builder.check_structure(wdir, original_species, protein) for protein in all_proteins]
        # missing_structures_final = [structure for structure in missing_structures if structure != None]

    # if user_input4 == "y" and missing_structures_final != []:
    #    print("...WARNING...: (FREEDA) I did not find clear structure prediction models for: %s" % missing_structures_final)
    #    print("...WARNING...: (FREEDA) I cannot overlay adaptive sites for these proteins\n")

    # ----------------------------------------#
    ######## RUN BLAST ########
    # ----------------------------------------#

    if user_input1 == "y":
        print("\n -> checking genome blast databases...")
        blast_path = tblastn.run_blast(wdir, original_species)
        if blast_path == None:
            print("\n...FATAL ERROR...: Blast database build failed for at least one genome"
                  "\n   Make sure you downloaded all genomes -> exiting the pipeline now...")
            return

    else:
        blast_path = wdir + "Blast_output/"

    # ADD A CHECK FOR BLAST OUTPUT IN CASE IT WAS TEMPERED WITH OR SPECIES WERE REMOVED

    # ----------------------------------------#
    ######## RUN EXON FINDING ########
    # ----------------------------------------#

    if user_input2 == "y":
        result_path = exon_extractor.analyse_blast_results(wdir, blast_path, original_species, int(t))

    # ----------------------------------------#
    ######## RUN PAML ########
    # ----------------------------------------#

    if user_input3 == "y" and user_input2 == "n":
        nr_of_tries = 1
        while nr_of_tries <= 3:
            user_input4 = input("----->  Indicate folder with results (no slashes or quotes): ")
            result_path = wdir + user_input4 + "/"
            if os.path.isdir(result_path) == False:
                print("\n(FREEDA) I could not find your results folder (%s/3)" % nr_of_tries)
                nr_of_tries += 1
            else:
                nr_of_tries = float("inf")
                result_path = wdir + user_input4 + "/"
                proteins, nr_of_species_total_dict, PAML_logfile_name, day = paml_launcher.analyse_final_cds(wdir,
                                                                                                             original_species,
                                                                                                             result_path)
                paml_visualizer.analyse_PAML_results(wdir, result_path,
                                                     proteins, nr_of_species_total_dict, original_species,
                                                     PAML_logfile_name, day)

    if user_input3 == "y" and user_input2 == "y":
        proteins, nr_of_species_total_dict, PAML_logfile_name, day = paml_launcher.analyse_final_cds(wdir,
                                                                                                     original_species,
                                                                                                     result_path)
        paml_visualizer.analyse_PAML_results(wdir, result_path,
                                             proteins, nr_of_species_total_dict, original_species, PAML_logfile_name,
                                             day)

    # ----------------------------------------#
    ######## RUN PyMOL ########
    # ----------------------------------------#

    # check_structure is obsolete cose not every model is usable -> refer to input_dictionary.txt

    if user_input1 == "n" and user_input2 == "n" and user_input3 == "n" and user_input4 == "y":
        nr_of_tries = 1
        while nr_of_tries <= 3:
            user_input4 = input("----->  Indicate folder with results (no slashes or quotes): ")
            result_path = wdir + user_input4 + "/"
            if os.path.isdir(result_path) == False:
                print("\n(FREEDA) I could not find your results folder (%s/3)" % nr_of_tries)
                nr_of_tries += 1
            else:
                # all_proteins = [protein.rstrip("\n") for protein in open(wdir + "proteins.txt", "r").readlines()]
                nr_of_tries = float("inf")
                result_path = wdir + user_input4 + "/"

                for protein, model_equal_input in input_dictionary.items():
                    # check if model seq and input seq match (from dict) and check if exactly one model exists
                    if bool(model_equal_input) == True and structure_builder.check_structure(wdir, original_species,
                                                                                             protein):
                        successful = structure_builder.run_pymol(wdir, original_species, result_path, protein,
                                                                 offset=None)
                        if not successful:
                            break
                    else:
                        print(
                            "\nPrediction model for : %s DOES NOT match input sequence -> cannot run PyMOL\n" % protein)

                # check if all models are present and run PyMOL
                # missing_structures_final = structure_builder.check_all_structures(wdir, original_species)
                # for protein in all_proteins:
                # check if model sequence equals the blasted sequence
                # model_equal_input = structure_builder.compare_model_with_input(wdir, original_species, protein, model_seq)
                # if protein not in missing_structures_final and model_equal_input == True:
                # structure_builder.run_pymol(wdir, original_species, result_path, protein, offset=None)

                # indicate that structure overlay is already present
                structure_overlay_present = True

    # run PyMOL together with the rest of the pipeline
    if user_input4 == "y" and not structure_overlay_present:
        for protein, model_equal_input in input_dictionary.items():
            # check if model sequence equals the blasted sequence
            # model_equal_input = structure_builder.compare_model_with_input(wdir, original_species, protein, model_seq)
            # check if model seq and input seq match (from dict) and check if exactly one model exists
            if bool(model_equal_input) == True and structure_builder.check_structure(wdir, original_species, protein):
                successful = structure_builder.run_pymol(wdir, original_species, result_path, protein, offset=None)
                if not successful:
                    break
            else:
                print("\nPrediction model for : %s DOES NOT match input sequence -> cannot run PyMOL\n" % protein)

    print("\nYou reached the end of FREEDA pipeline.")


if __name__ == '__freeda_pipeline__':
    import argparse

    parser = argparse.ArgumentParser(description="Run FREEDA pipeline")
    parser.add_argument("-os", "--original_species",
                        help="specify reference organism (ex. Mm for Mus musculus)", type=str,
                        required=False)
    parser.add_argument("-t", "--blast_threshold",
                        help="specify percentage identity threshold for blast (ex. 30; default)", type=int,
                        required=True)

    args = parser.parse_args()
    freeda_pipeline(original_species=args.original_species, t=args.blast_threshold)

"""

 
                # get structure model using AlphaFold url request
                prediction_url = input_extractor.get_prediction(wdir, original_species, protein)
                if prediction_url == None:
                    print("AlphaFold prediction not available for: %s\n" % protein)
                    model_equal_input = False
                    pass
                elif prediction_url == True:
                    print("Structure prediction model for: %s already exists\n" % protein)
                    # check if model sequence equals the blasted sequence
                    model_equal_input = structure_builder.compare_model_with_input(wdir, original_species, protein)
                    pass
                else:
                    print("\n(FREEDA) Please input structure prediction for protein: %s\n(copy the following url into your browser " \
                      "-> click PDB file -> save in ../Data/Structures/%s)\n\n " \
                         "%s\n\n ...WARNING... Verify protein identity (if incorrect find model in AlphaFold browser)" 
                                                     % (protein, protein + "_" + original_species, prediction_url))
                    input("\n(FREEDA) When done press ENTER\n")
                    # check if model sequence equals the blasted sequence
                    model_equal_input = structure_builder.compare_model_with_input(wdir, original_species, protein)


"""
