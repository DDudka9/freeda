#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:36:16 2021

@author: damian
Main module of the freeda package. Takes user input and performs automatic input extraction, tblastn, exon finding
and molecular evolution analysis (PAML) followed by overlay of putative adaptive sites onto 3D structure (PyMOL).

"""


# TODO
#       0) cite Wang and Han 2021 J Virol for Primates, Carnivora that span similar phylogeny
#       "Pervasive Positive Selection on Virus Receptors Driven by Host-Virus Conflicts in Mammals"
#       0) There might be an issue with RETRO calling -> Only 5-prime RETRO called in Mad1l1 Mi (the one that didnt fail)
#       0) Whn only F61 scores (NUMA1) PAML log has annotated uniprot like sites in F3X4 too -> need to fix it
#               -> PAML graph also shows sites in NUMA1 F3X4 model -> DONE (I think ->test on NUMA1)
#                   -> it works but structure overlays only F61 sites that also score in F3X4
#                                                                       (even if this model is rejected)
#               -> it turns out that actually freeda picks the first codon model if consensus dict is not neded (fix it)
#       0) Test current folder scheme using command line
#       0) When F61 scores only, then structure has these residues -> its ok, command line can run only F61
#       0) Cenpk using Rn as reference, Ha genome exon 8 is eliminated but introns are 0.74 and 1.0
#                   -> the 1.0 intron is actually completely missing -> no mismatches are treated as full alignment
#       0) You need to introduce info about the ref exon length and pass it to single exon checkpoint to avoid
#               the code recognising intronic dashes as indels -> one of the dicts "ref_exons" should have it
#                   -> I added the exon length variable to single exon mapping checkpoint, penalized gaps
#                           by giving 0.25 point instead of 1 -> test on CD55
#                   -> I decided to go back to not penalizing indels because that leads to rejection of exons e.g.
#                           Cenp-a ends up with only 7 species becuase of the exon 1 indels (which are true)
#       0) CD55 in primates An species has clear misalignment in exon 4 and 5 -> passes single exon check though (0.69 and 0.74)
#                   -> I started punishing indels -> test on Izumo1, Izumo 3 and Haus8 -> decided against it finally?
#       0) Possible issue with length of gene -> how TRIM5a human gene is soo long? Trim-214 in ensembl is
#                   only 20k while I get over 200k linear seq -> check input from pyensembl
#                               -> seems that its pyensembl doing; gives bed coordinates for gene that long -< no clue why
#       0) Consider papers when writing: "The effects of alignment error and alignment filtering on the sitewise detection of positive selection"
#                                       : "A beginners guide to estimating the non-synonymous to
#                                                 synonymous rate ratio of all protein-coding genes in a genome"
#       0) Program to make "denovo codon-based nucleotide alignment" when protein-coding sequences are unreliable or absent":
#                   -> MACSE "Ranwez et al., 2011 PLoS One"
#       0) Jeffares et al., 2015 Methods Mol Biol A beginners guide to estimating the non-synonymous to
#               synonymous rate ratio of all protein-coding genes in a genome
#                           -> indentity of 70% produces reliable alignments
#       0) Lccr37a -> leucine rich repetitive protein -> example of protein FREEDA cannot analyse
#       0) Get rid of the "another intron interpheres with introny check" for clarity
#       0) Consider showing only 0.95 adaptive sites on structure if more than 20 sites
#                       or only 0.99 when more than 50 sites -> not sure
#       0) Either disable the interpro domains or restrict them even more (too overlapping) -> not sure
#       0) Change pr >= 0.9 residues to black if 0 -> DONE
#       0) Make a file from the list of genomes -> will allow users to add more genomes are more is being sequenced
#                   -> DONE
#       0) Similar to Cenpn, terf2 rat is only 92% similar to ensembl but 100% identical to uniprot predited sequence
#                   -> disable the check in the final freeda verison
#       0) General issue -> Traceback are not logged in, in case of crashing the GUI will not stop
#                               -> confusing for user
#                           SOLUTION : find a way to crash the GUI when exception occurs that FREEDA does not handle
#                           -> tried to add it as a try/except statement -> doesnt work
#       0) Figures:
#               1) Accuracy - Show both ways (rat and mouse)
#                           - Show mafft vs muscle
#                           - Show dealing with duplications (Pot1a and Pot1b?)
#                           - Show dealing with tandem sequences (Hoxd9, Hoxd10)
#                                   -> warning that it cannot distinguish tandem duplications
#                                               that are too recent, especially in lower quality genomes (Mug1, Mug2)
#                           - Show dealing with large introns (Cenpp)
#                           - Show dealing with microexons (Cenpx)
#                           - Show dealing with premature STOP (Haus...)
#                           - Show dealing with uncalled bases (Mug1 Caroli contig FMAL02029158.1__rev)
#       0) Use debugger to go through the matches_processor module and make sure docstrings are correct
#       0) Double check all the assemblies (I swapped white faced saki in order)
#       2) Talk to Mike about who to ask concering licenses
#       3) Fix nomenclature -> "intronic" -> DONE
#                           -> "gene/gene_name" -> gene_name -> DONE
#                           -> "None" is not the best info for M2a etc tests in PAML excel file -> 0 (zero)
#                                   -> DONE (testing...)
#                           -> eliminate aligner option (at the end when figure is done)
#                           -> fill out missing docstrings
#      12) Use Hyland et al. 2021 Gen Biol Evol for testing -> TRIP gene in mammals

from freeda import input_extractor
from freeda import folder_generator
from freeda import tblastn
from freeda import exon_extractor
from freeda import paml_launcher
from freeda import paml_visualizer
from freeda import structure_builder
from freeda import genomes_preprocessing
from freeda import TextHandler
from freeda import pyinstaller_compatibility
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from PIL import ImageTk, Image
import tkinter.scrolledtext as ScrolledText
import os
import shutil
import re
import logging
import threading
from sys import platform


def thread_freeda():
    """Runs freeda main function in another thread"""

    if check_input():

        freeda_thread = threading.Thread(target=freeda_pipeline)
        # makes a daemon thread allowing it to be killed anytime (though not advised generally)
        freeda_thread.daemon = True
        freeda_thread.start()
        block_user_entries()
        # erase previous results
        clear_results()

    else:
        logging.info("\n--------------- TRY AGAIN :D -------------------\n")
        # unblock user input entries
        ublock_user_entries()


def abort_freeda():
    """Closes the GUI and aborts FREEDA pipeline"""

    root.quit()
    root.destroy()


def raise_logger():
    """Raises a logger for the GUI"""

    global logging_window

    # LOGGER
    logging_label = ttk.Label(logging_frame, text="Events window (logged to 'FREEDA*.log' and 'PAML*.log')")
    logging_label.grid(column=0, row=0, columnspan=4, sticky=(W))
    logging_window = ScrolledText.ScrolledText(logging_frame, state="disabled", wrap="none")
    logging_window.grid(column=0, row=1, columnspan=7, sticky=(N, W, E, S))
    logging_window.configure(font='TkFixedFont')

    # create handlers of the logging window
    text_handler = TextHandler.TextHandler(logging_window)
    # Logging configuration
    logging.basicConfig(format="%(message)s")  # level=logging.INFO,
    logger = logging.getLogger()
    logger.addHandler(text_handler)


def check_input():
    """Checks for essential user input"""

    ready = True

    if not wdirectory.get():
        logging.info("\n...FATAL_ERROR... : Choose working directory")
        ready = False

    if not clade.get():
        logging.info("\n...FATAL_ERROR... : Choose clade")
        ready = False

    all_genes = [gene_name1.get(), gene_name2.get(), gene_name3.get(), gene_name4.get(), gene_name5.get()]
    if not any(all_genes):
        logging.info("\n...FATAL_ERROR... : Choose at least one gene")
        ready = False

    # check if there are duplications in user input (but not empty entries)
    if len(all_genes) != len(set(all_genes)):
        if gene_name1.get() == gene_name2.get() and gene_name1.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False
        elif gene_name1.get() == gene_name3.get() and gene_name3.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False
        elif gene_name1.get() == gene_name4.get() and gene_name1.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False
        elif gene_name1.get() == gene_name5.get() and gene_name1.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False
        elif gene_name2.get() == gene_name3.get() and gene_name2.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False
        elif gene_name2.get() == gene_name4.get() and gene_name2.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False
        elif gene_name2.get() == gene_name5.get() and gene_name2.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False
        elif gene_name3.get() == gene_name4.get() and gene_name3.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False
        elif gene_name3.get() == gene_name5.get() and gene_name3.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False
        elif gene_name4.get() == gene_name5.get() and gene_name4.get() != "":
            logging.info("\n...FATAL_ERROR... : Choose different gene names")
            ready = False

    return ready


def block_user_entries():
    """Blocks the user input entries"""

    # disable user input entries
    rodents.configure(state="disabled")
    primates.configure(state="disabled")
    carnivores.configure(state="disabled")
    birds.configure(state="disabled")

    codon_freq_F3X4.configure(state="disabled")
    codon_freq_F61.configure(state="disabled")

    name1.configure(state="disabled")
    name2.configure(state="disabled")
    name3.configure(state="disabled")
    name4.configure(state="disabled")
    name5.configure(state="disabled")

    dup1_button.configure(state="disabled")
    dup2_button.configure(state="disabled")
    dup3_button.configure(state="disabled")
    dup4_button.configure(state="disabled")
    dup5_button.configure(state="disabled")

    site11_start.configure(state="disabled")
    site11_end.configure(state="disabled")
    site11_label.configure(state="disabled")
    site12_start.configure(state="disabled")
    site12_end.configure(state="disabled")
    site12_label.configure(state="disabled")
    site13_start.configure(state="disabled")
    site13_end.configure(state="disabled")
    site13_label.configure(state="disabled")

    site21_start.configure(state="disabled")
    site21_end.configure(state="disabled")
    site21_label.configure(state="disabled")
    site22_start.configure(state="disabled")
    site22_end.configure(state="disabled")
    site22_label.configure(state="disabled")
    site23_start.configure(state="disabled")
    site23_end.configure(state="disabled")
    site23_label.configure(state="disabled")

    site31_start.configure(state="disabled")
    site31_end.configure(state="disabled")
    site31_label.configure(state="disabled")
    site32_start.configure(state="disabled")
    site32_end.configure(state="disabled")
    site32_label.configure(state="disabled")
    site33_start.configure(state="disabled")
    site33_end.configure(state="disabled")
    site33_label.configure(state="disabled")

    site41_start.configure(state="disabled")
    site41_end.configure(state="disabled")
    site41_label.configure(state="disabled")
    site42_start.configure(state="disabled")
    site42_end.configure(state="disabled")
    site42_label.configure(state="disabled")
    site43_start.configure(state="disabled")
    site43_end.configure(state="disabled")
    site43_label.configure(state="disabled")

    site51_start.configure(state="disabled")
    site51_end.configure(state="disabled")
    site51_label.configure(state="disabled")
    site52_start.configure(state="disabled")
    site52_end.configure(state="disabled")
    site52_label.configure(state="disabled")
    site53_start.configure(state="disabled")
    site53_end.configure(state="disabled")
    site53_label.configure(state="disabled")

    exclude_species_entry.configure(state="disabled")

    wdir_button.configure(state="disabled")
    wdir_entry.configure(state="disabled")


def ublock_user_entries():
    """Unblocks the disabled entries allowing running FREEDA again"""

    # disable user input entries
    rodents.configure(state="normal")
    primates.configure(state="normal")
    carnivores.configure(state="normal")
    birds.configure(state="normal")

    codon_freq_F3X4.configure(state="normal")
    codon_freq_F61.configure(state="normal")

    name1.configure(state="normal")
    name2.configure(state="normal")
    name3.configure(state="normal")
    name4.configure(state="normal")
    name5.configure(state="normal")

    dup1_button.configure(state="normal")
    dup2_button.configure(state="normal")
    dup3_button.configure(state="normal")
    dup4_button.configure(state="normal")
    dup5_button.configure(state="normal")

    site11_start.configure(state="normal")
    site11_end.configure(state="normal")
    site11_label.configure(state="normal")
    site12_start.configure(state="normal")
    site12_end.configure(state="normal")
    site12_label.configure(state="normal")
    site13_start.configure(state="normal")
    site13_end.configure(state="normal")
    site13_label.configure(state="normal")

    site21_start.configure(state="normal")
    site21_end.configure(state="normal")
    site21_label.configure(state="normal")
    site22_start.configure(state="normal")
    site22_end.configure(state="normal")
    site22_label.configure(state="normal")
    site23_start.configure(state="normal")
    site23_end.configure(state="normal")
    site23_label.configure(state="normal")

    site31_start.configure(state="normal")
    site31_end.configure(state="normal")
    site31_label.configure(state="normal")
    site32_start.configure(state="normal")
    site32_end.configure(state="normal")
    site32_label.configure(state="normal")
    site33_start.configure(state="normal")
    site33_end.configure(state="normal")
    site33_label.configure(state="normal")

    site41_start.configure(state="normal")
    site41_end.configure(state="normal")
    site41_label.configure(state="normal")
    site42_start.configure(state="normal")
    site42_end.configure(state="normal")
    site42_label.configure(state="normal")
    site43_start.configure(state="normal")
    site43_end.configure(state="normal")
    site43_label.configure(state="normal")

    site51_start.configure(state="normal")
    site51_end.configure(state="normal")
    site51_label.configure(state="normal")
    site52_start.configure(state="normal")
    site52_end.configure(state="normal")
    site52_label.configure(state="normal")
    site53_start.configure(state="normal")
    site53_end.configure(state="normal")
    site53_label.configure(state="normal")

    exclude_species_entry.configure(state="normal")

    wdir_button.configure(state="normal")
    wdir_entry.configure(state="normal")


def clear_results():
    """Erases the results from the previous run"""

    g1_results_entry.configure(state="normal")
    g1_results_entry.delete(0, 'end')
    g1_results_entry.configure(state="disabled")
    g1_pos_sel_entry.configure(state="normal")
    g1_pos_sel_entry.delete(0, 'end')
    g1_pos_sel_entry.configure(state="disabled")
    g1_lrt_entry.configure(state="normal")
    g1_lrt_entry.delete(0, 'end')
    g1_lrt_entry.configure(state="disabled")
    g1_pvalue_entry.configure(state="normal")
    g1_pvalue_entry.delete(0, 'end')
    g1_pvalue_entry.configure(state="disabled")
    g1_coverage_entry.configure(state="normal")
    g1_coverage_entry.delete(0, 'end')
    g1_coverage_entry.configure(state="disabled")
    g1_species_entry.configure(state="normal")
    g1_species_entry.delete(0, 'end')
    g1_species_entry.configure(state="disabled")
    g1_adapt_less_entry.configure(state="normal")
    g1_adapt_less_entry.delete(0, 'end')
    g1_adapt_less_entry.configure(state="disabled")
    g1_adapt_more_entry.configure(state="normal")
    g1_adapt_more_entry.delete(0, 'end')
    g1_adapt_more_entry.configure(state="disabled")

    g2_results_entry.configure(state="normal")
    g2_results_entry.delete(0, 'end')
    g2_results_entry.configure(state="disabled")
    g2_pos_sel_entry.configure(state="normal")
    g2_pos_sel_entry.delete(0, 'end')
    g2_pos_sel_entry.configure(state="disabled")
    g2_lrt_entry.configure(state="normal")
    g2_lrt_entry.delete(0, 'end')
    g2_lrt_entry.configure(state="disabled")
    g2_pvalue_entry.configure(state="normal")
    g2_pvalue_entry.delete(0, 'end')
    g2_pvalue_entry.configure(state="disabled")
    g2_coverage_entry.configure(state="normal")
    g2_coverage_entry.delete(0, 'end')
    g2_coverage_entry.configure(state="disabled")
    g2_species_entry.configure(state="normal")
    g2_species_entry.delete(0, 'end')
    g2_species_entry.configure(state="disabled")
    g2_adapt_less_entry.configure(state="normal")
    g2_adapt_less_entry.delete(0, 'end')
    g2_adapt_less_entry.configure(state="disabled")
    g2_adapt_more_entry.configure(state="normal")
    g2_adapt_more_entry.delete(0, 'end')
    g2_adapt_more_entry.configure(state="disabled")

    g3_results_entry.configure(state="normal")
    g3_results_entry.delete(0, 'end')
    g3_results_entry.configure(state="disabled")
    g3_pos_sel_entry.configure(state="normal")
    g3_pos_sel_entry.delete(0, 'end')
    g3_pos_sel_entry.configure(state="disabled")
    g3_lrt_entry.configure(state="normal")
    g3_lrt_entry.delete(0, 'end')
    g3_lrt_entry.configure(state="disabled")
    g3_pvalue_entry.configure(state="normal")
    g3_pvalue_entry.delete(0, 'end')
    g3_pvalue_entry.configure(state="disabled")
    g3_coverage_entry.configure(state="normal")
    g3_coverage_entry.delete(0, 'end')
    g3_coverage_entry.configure(state="disabled")
    g3_species_entry.configure(state="normal")
    g3_species_entry.delete(0, 'end')
    g3_species_entry.configure(state="disabled")
    g3_adapt_less_entry.configure(state="normal")
    g3_adapt_less_entry.delete(0, 'end')
    g3_adapt_less_entry.configure(state="disabled")
    g3_adapt_more_entry.configure(state="normal")
    g3_adapt_more_entry.delete(0, 'end')
    g3_adapt_more_entry.configure(state="disabled")

    g4_results_entry.configure(state="normal")
    g4_results_entry.delete(0, 'end')
    g4_results_entry.configure(state="disabled")
    g4_pos_sel_entry.configure(state="normal")
    g4_pos_sel_entry.delete(0, 'end')
    g4_pos_sel_entry.configure(state="disabled")
    g4_lrt_entry.configure(state="normal")
    g4_lrt_entry.delete(0, 'end')
    g4_lrt_entry.configure(state="disabled")
    g4_pvalue_entry.configure(state="normal")
    g4_pvalue_entry.delete(0, 'end')
    g4_pvalue_entry.configure(state="disabled")
    g4_coverage_entry.configure(state="normal")
    g4_coverage_entry.delete(0, 'end')
    g4_coverage_entry.configure(state="disabled")
    g4_species_entry.configure(state="normal")
    g4_species_entry.delete(0, 'end')
    g4_species_entry.configure(state="disabled")
    g4_adapt_less_entry.configure(state="normal")
    g4_adapt_less_entry.delete(0, 'end')
    g4_adapt_less_entry.configure(state="disabled")
    g4_adapt_more_entry.configure(state="normal")
    g4_adapt_more_entry.delete(0, 'end')
    g4_adapt_more_entry.configure(state="disabled")

    g5_results_entry.configure(state="normal")
    g5_results_entry.delete(0, 'end')
    g5_results_entry.configure(state="disabled")
    g5_pos_sel_entry.configure(state="normal")
    g5_pos_sel_entry.delete(0, 'end')
    g5_pos_sel_entry.configure(state="disabled")
    g5_lrt_entry.configure(state="normal")
    g5_lrt_entry.delete(0, 'end')
    g5_lrt_entry.configure(state="disabled")
    g5_pvalue_entry.configure(state="normal")
    g5_pvalue_entry.delete(0, 'end')
    g5_pvalue_entry.configure(state="disabled")
    g5_coverage_entry.configure(state="normal")
    g5_coverage_entry.delete(0, 'end')
    g5_coverage_entry.configure(state="disabled")
    g5_species_entry.configure(state="normal")
    g5_species_entry.delete(0, 'end')
    g5_species_entry.configure(state="disabled")
    g5_adapt_less_entry.configure(state="normal")
    g5_adapt_less_entry.delete(0, 'end')
    g5_adapt_less_entry.configure(state="disabled")
    g5_adapt_more_entry.configure(state="normal")
    g5_adapt_more_entry.delete(0, 'end')
    g5_adapt_more_entry.configure(state="disabled")


def freeda_pipeline():
    """Main function running all freeda pipeline"""

    # main function is enclosed in a try/except to avoid frozen GUI
    try:
        # use of the pipeline from GUI
        gui = True

        # deactivate the analyze button
        analyze_button["state"] = "disable"

        # get user input
        wdir = wdirectory.get() + "/"
        ref_species = clade.get()
        codon_frequencies = codon_freq.get()
        t = 60
        all_genes_dict = {gene_name1.get(): [dup1_var.get(),
                                        [site11_label.get(), site11_start.get(), site11_end.get()],
                                        [site12_label.get(), site12_start.get(), site12_end.get()],
                                        [site13_label.get(), site13_start.get(), site13_end.get()]],
                             gene_name2.get(): [dup2_var.get(),
                                        [site21_label.get(), site21_start.get(), site21_end.get()],
                                        [site22_label.get(), site22_start.get(), site22_end.get()],
                                        [site23_label.get(), site23_start.get(), site23_end.get()]],
                             gene_name3.get(): [dup3_var.get(),
                                        [site31_label.get(), site31_start.get(), site31_end.get()],
                                        [site32_label.get(), site32_start.get(), site32_end.get()],
                                        [site33_label.get(), site33_start.get(), site33_end.get()]],
                             gene_name4.get(): [dup4_var.get(),
                                        [site41_label.get(), site41_start.get(), site41_end.get()],
                                        [site42_label.get(), site42_start.get(), site42_end.get()],
                                        [site43_label.get(), site43_start.get(), site43_end.get()]],
                             gene_name5.get(): [dup5_var.get(),
                                        [site51_label.get(), site51_start.get(), site51_end.get()],
                                        [site52_label.get(), site52_start.get(), site52_end.get()],
                                        [site53_label.get(), site53_start.get(), site53_end.get()]]}
        # get a list of genes
        all_genes = [gene for gene in all_genes_dict if gene != ""]

        excluded_species = exclude_species_var.get()
        excluded_species = excluded_species.split(" ")
        final_excluded_species = {}
        available_species = genomes_preprocessing.get_available_species(ref_species)
        for ex_species in excluded_species:
            # check if users input matches available species
            if ex_species in available_species:
                final_excluded_species[ex_species] = available_species[ex_species]  # add genome as value to species key

        global logging_window

        # LOGGER
        logging_label = ttk.Label(logging_frame, text="Events window (logged to 'FREEDA*.log' and 'PAML*.log')")
        logging_label.grid(column=0, row=0, columnspan=4, sticky=(W))
        logging_window = ScrolledText.ScrolledText(logging_frame, state="disabled", wrap="none")
        logging_window.grid(column=0, row=1, columnspan=7, sticky=(N, W, E, S))
        logging_window.configure(font='TkFixedFont')

        # create handlers of the logging window
        text_handler = TextHandler.TextHandler(logging_window)
        # Logging configuration
        logging.basicConfig(format="%(message)s")  # level=logging.INFO,
        logger = logging.getLogger()
        logger.addHandler(text_handler)

        # change working directory to user indicated
        os.chdir(wdir)

        # ----------------------------------------#
        ######## GET USER INPUT ########
        # ----------------------------------------#

        # generate basic folders for input if not present
        folder_generator.generate_basic_folders(wdir)

        # get settings
        aligner = "mafft"

        # get all species and genome names
        all_genomes = [genome[1] for genome in genomes_preprocessing.get_names(wdir, ref_species)]

        # ----------------------------------------#
        ######## INSTALL PYMOL IF NEEDED ########
        # ----------------------------------------#

        if not shutil.which("pymol"):
            if platform == "linux" or platform == "linux2":
                logging.info("\nPyMOL not found in the PATH. Installing PyMOL to the Data folder now.")
                structure_builder.install_pymol_linux(wdir)
            else:
                print("PyMOL")
            return

        # ----------------------------------------#
        ######## GET ALL INPUT DATA  ########
        # ----------------------------------------#

        # generate a reference Genome object
        ref_genome_present, ensembl, ref_species, ref_genomes_path, ref_genome_contigs_dict, \
            biotype, all_genes_ensembl = input_extractor.generate_ref_genome_object(wdir, ref_species)

        if not input_extractor.validate_gene_names(all_genes, all_genes_ensembl):
            ublock_user_entries()
            return

        # stop pipeline if the reference genome is absent
        if not ref_genome_present:
            message = "\n...FATAL ERROR... : There is no reference genome detected -> exiting the pipeline now...\n"
            logging.info(message)
            ublock_user_entries()
            return

        # get names of genes

        for gene in all_genes:

            message = "\n----------- * %s * -----------" % gene
            logging.info(message)
            # get structure prediction model from AlphaFold
            possible_uniprot_ids = input_extractor.get_uniprot_id(ref_species, gene)
            model_seq, uniprot_id = input_extractor.fetch_structure_prediction(wdir, ref_species,
                                                                                   gene, possible_uniprot_ids)
            # get sequence input from ensembl
            input_correct, model_matches_input, microexon_present, microexons = input_extractor.extract_input(
                wdir, ref_species, ref_genomes_path,
                ref_genome_contigs_dict, ensembl, biotype,
                gene, model_seq, uniprot_id
            )

            if input_correct:
                message = "\nInput data have been generated for gene: %s\n\n" % gene
                logging.info(message)

            if not input_correct:
                message = "\n...FATAL ERROR... : Input data generation FAILED for gene: %s " \
                          "- please remove from analysis -> exiting the pipeline now...\n" % gene
                logging.info(message)
                ublock_user_entries()
                return

            if not model_matches_input:
                message = "\n...WARNING... : No matching structure prediction model is available for : %s " \
                      "-> cannot overlay FREEDA results onto a 3D structure\n" % gene
                logging.info(message)

        # ----------------------------------------#
        ######## RUN BLAST ########
        # ----------------------------------------#

        print("Checking genome blast databases...")
        tblastn.run_blast(wdir, ref_species, all_genes, final_excluded_species)

        # ----------------------------------------#
        ######## RUN EXON FINDING ########
        # ----------------------------------------#

        if exon_extractor.check_blast_output(wdir + "Blast_output/", t, all_genes):
            result_path = exon_extractor.analyze_blast_results(wdir, wdir + "Blast_output/",
                                                               ref_species, int(t), all_genes, all_genomes, aligner,
                                                               final_excluded_species, gui, logging_window,
                                                               all_genes_dict)
            # set a StringVar for GUI
            result_path_var.set(result_path)

        else:
            message = "\n     ...FATAL ERROR... : Genome of at least one species contains " \
                  "no matches above the identity threshold used : %s \n" \
                      "-> use a lower one or exclude from analysis" \
                  "-> exiting the pipeline now..." % t
            logging.info(message)
            ublock_user_entries()
            return

        # ----------------------------------------#
        ######## RUN PAML and PyMOL ########
        # ----------------------------------------#

        # run PAML
        nr_of_species_total_dict, PAML_logfile_name, day, failed_paml, \
                genes_under_pos_sel = paml_launcher.analyze_final_cds(wdir, ref_species, result_path,
                                                                      all_genes, aligner, codon_frequencies,
                                                                      gui, logging_window)

        # visualize PAML result
        final_PAML_log_dict = paml_visualizer.analyze_PAML_results(wdir, result_path, all_genes,
                                                                   nr_of_species_total_dict, ref_species,
                                                                   PAML_logfile_name, day, genes_under_pos_sel,
                                                                   failed_paml, codon_frequencies, gui)
        # in case PAML failed and dict wasnt created
        if final_PAML_log_dict:
            get_results(final_PAML_log_dict)
        else:
            ublock_user_entries()
            return

        # run PyMOL
        for gene in all_genes:

            # do not allow further analysis of failed paml runs
            if gene in failed_paml:
                continue

            # check if model seq and input seq match and check if exactly one model exists
            elif structure_builder.check_structure(wdir, ref_species, gene):
                successful = structure_builder.run_pymol(wdir, ref_species, result_path, gene, genes_under_pos_sel,
                                                                     all_genes_dict)
                if not successful:
                    message = "\nThe structure for : %s was not built successfully." % gene
                    logging.info(message)
                    continue
            else:
                message = "\nPrediction model for : %s DOES NOT match input sequence" \
                         "-> cannot overlay FREEDA results onto a 3D structure\n" % gene
                logging.info(message)

        logging.info("\nYou reached the end of FREEDA pipeline.")

        # allow user to run freeda again
        ublock_user_entries()

        # reset the logger to avoid writing into the previous log file
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

    # report exceptions
    except Exception:
        logging.exception("\n\n...FATAL ERROR... : Something went wrong (see below). "
                          "Send screen shot to : Damian Dudka -> damiandudka0@gmail.com\n\n")


def check_gene_name(gene_name, op):
    """Checks if user provided a valid gene name"""
    error_message1.set("")
    # accept only entry starting with one capital letter followed by small or big letters or numbers
    valid = re.match(r"^[A-Z]{1}([A-Za-z0-9]+$)", gene_name) is not None
    # button can be clicked only if gene names are valid
    analyze_button.state(["!disabled"] if valid else ["disabled"])
    # keystroke validation
    if op == "key":
        ok_so_far = re.match(r"^(?![\s\S])|[A-Za-z0-9]+$", gene_name) is not None  # ^(?![\s\S]) -> completely empty
        if not ok_so_far:
            error_message1.set(message1)
        return ok_so_far
    elif op == "focusout":
        if not valid:
            error_message1.set(message1)
    return valid


def check_functional_residues(residue, op):
    """Checks if user provided a valid residue number"""
    error_message2.set("")
    # accept only entry starting with non-0 number followed by max 3 numbers
    valid = re.match(r"^[1-9]{1}([0-9]{0,3}$)", residue) is not None  # max 4 digits, first one not a zero
    # keystroke validation
    if op == "key":
        ok_so_far = re.match(r"^(?![\s\S])|[0-9]+$", residue) is not None  # ^(?![\s\S]) -> completely empty
        if not ok_so_far:
            error_message2.set(message2)
        return ok_so_far
    elif op == "focusout":
        if not valid:
            error_message2.set(message2)
    return valid


def check_label(label, op):
    """Checks if user provided a valid label name"""
    # accept only entry composed of letters and numbers
    valid = re.match(r"^[A-Za-z0-9]+$", label) is not None
    # keystroke validation
    if op == "key":
        ok_so_far = re.match(r"^(?![\s\S])|[\w\s]+$", label) is not None  # ^(?![\s\S]) -> completely empty
        return ok_so_far

    return valid


def check_excluded_species(label, op):
    """Checks if user provided a valid species name"""
    # accept only entry composed of letters and numbers
    valid = re.match(r"^[A-Z]{1}[a-z]{1}$", label) is not None
    # keystroke validation
    if op == "key":
        ok_so_far = re.match(r"^(?![\s\S])|[\w\s]+$", label) is not None  # ^(?![\s\S]) -> completely empty
        return ok_so_far

    return valid


def get_wdir():
    """Asks for a working directory."""
    directory = filedialog.askdirectory(initialdir=os.getcwd(), title="Select 'Data' folder")
    wdirectory.set(directory)


def get_results(final_PAML_log_dict):
    """Fetches key PAML results for each gene"""

    genes = final_PAML_log_dict["Gene name"]
    lrts = final_PAML_log_dict["M8 vs M7 (LRT)"]
    pvalues_m7m8 = final_PAML_log_dict["M8 vs M7 (p-value)"]
    pvalues_m1m2 = final_PAML_log_dict["M2a vs M1a (p-value)"]
    coverage = final_PAML_log_dict["CDS Coverage"]
    species = final_PAML_log_dict["Nr of species analyzed"]
    adapt_less = final_PAML_log_dict["Sites with pr < 0.90"]
    adapt_more = final_PAML_log_dict["Sites with pr >= 0.90"]

    if gene_name1.get():

        g1_results_var.set(genes.pop(0))
        lrt = lrts.pop(0)
        if lrt == "NA":
            g1_lrt_var.set("NA")
        elif lrt == "(0)":
            g1_lrt_var.set("(0)")
        else:
            g1_lrt_var.set(round(float(lrt), ndigits=4))

        pvalue_m7m8 = pvalues_m7m8.pop(0)
        pvalue_m1m2 = pvalues_m1m2.pop(0)
        if pvalue_m7m8 == "NA":
            g1_pvalue_var.set("NA")
            g1_pos_sel_var.set("NA")
            g1_pos_sel_entry.config(foreground="black")
        else:
            pvalue_m7m8 = float(pvalue_m7m8)
            pvalue_m1m2 = float(pvalue_m1m2)
            if pvalue_m7m8 < 0.05:
                if pvalue_m1m2 < 0.05:
                    g1_pos_sel_var.set("YES")
                else:
                    g1_pos_sel_var.set("(YES)")
                g1_pos_sel_entry.config(foreground="magenta")
            if pvalue_m7m8 <= 0.001:
                g1_pvalue_var.set("<0.001")
            if pvalue_m7m8 > 0.001:
                g1_pvalue_var.set(round(float(pvalue_m7m8), ndigits=4))
            if pvalue_m7m8 >= 0.05:
                g1_pos_sel_var.set("NO")
                g1_pos_sel_entry.config(foreground="black")

        g1_coverage_var.set(coverage.pop(0))
        g1_species_var.set(species.pop(0))

        adapt_less_number = adapt_less.pop(0).split(" ")[0]
        if adapt_less_number == "Not":
            adapt_less_number = 0
        g1_adapt_less_var.set(adapt_less_number)

        adapt_more_number = adapt_more.pop(0).split(" ")[0]
        if adapt_more_number == "Not":
            adapt_more_number = 0
            g1_adapt_more_entry.config(foreground="black")
        elif adapt_more_number == "NA":
            g1_adapt_more_entry.config(foreground="black")
        else:
            g1_adapt_more_entry.config(foreground="magenta")
        g1_adapt_more_var.set(adapt_more_number)


    if gene_name2.get():
        g2_results_var.set(genes.pop(0))
        lrt = lrts.pop(0)
        if lrt == "NA":
            g2_lrt_var.set("NA")
        elif lrt == "(0)":
            g2_lrt_var.set("(0)")
        else:
            g2_lrt_var.set(round(float(lrt), ndigits=4))

        pvalue_m7m8 = pvalues_m7m8.pop(0)
        pvalue_m1m2 = pvalues_m1m2.pop(0)
        if pvalue_m7m8 == "NA":
            g2_pvalue_var.set("NA")
            g2_pos_sel_var.set("NA")
            g2_pos_sel_entry.config(foreground="black")
        else:
            pvalue_m7m8 = float(pvalue_m7m8)
            pvalue_m1m2 = float(pvalue_m1m2)
            if pvalue_m7m8 < 0.05:
                if pvalue_m1m2 < 0.05:
                    g2_pos_sel_var.set("YES")
                else:
                    g2_pos_sel_var.set("(YES)")
                g2_pos_sel_entry.config(foreground="magenta")
            if pvalue_m7m8 <= 0.001:
                g2_pvalue_var.set("<0.001")
            if pvalue_m7m8 > 0.001:
                g2_pvalue_var.set(round(float(pvalue_m7m8), ndigits=4))
            if pvalue_m7m8 >= 0.05:
                g2_pos_sel_var.set("NO")
                g2_pos_sel_entry.config(foreground="black")

        g2_coverage_var.set(coverage.pop(0))
        g2_species_var.set(species.pop(0))

        adapt_less_number = adapt_less.pop(0).split(" ")[0]
        if adapt_less_number == "Not":
            adapt_less_number = 0
        g2_adapt_less_var.set(adapt_less_number)

        adapt_more_number = adapt_more.pop(0).split(" ")[0]
        if adapt_more_number == "Not":
            adapt_more_number = 0
            g2_adapt_more_entry.config(foreground="black")
        elif adapt_more_number == "NA":
            g2_adapt_more_entry.config(foreground="black")
        else:
            g2_adapt_more_entry.config(foreground="magenta")
        g2_adapt_more_var.set(adapt_more_number)

    if gene_name3.get():
        g3_results_var.set(genes.pop(0))
        lrt = lrts.pop(0)
        if lrt == "NA":
            g3_lrt_var.set("NA")
        elif lrt == "(0)":
            g3_lrt_var.set("(0)")
        else:
            g3_lrt_var.set(round(float(lrt), ndigits=4))

        pvalue_m7m8 = pvalues_m7m8.pop(0)
        pvalue_m1m2 = pvalues_m1m2.pop(0)
        if pvalue_m7m8 == "NA":
            g3_pvalue_var.set("NA")
            g3_pos_sel_var.set("NA")
            g3_pos_sel_entry.config(foreground="black")
        else:
            pvalue_m7m8 = float(pvalue_m7m8)
            pvalue_m1m2 = float(pvalue_m1m2)
            if pvalue_m7m8 < 0.05:
                if pvalue_m1m2 < 0.05:
                    g3_pos_sel_var.set("YES")
                else:
                    g3_pos_sel_var.set("(YES)")
                g3_pos_sel_entry.config(foreground="magenta")
            if pvalue_m7m8 <= 0.001:
                g3_pvalue_var.set("<0.001")
            if pvalue_m7m8 > 0.001:
                g3_pvalue_var.set(round(float(pvalue_m7m8), ndigits=4))
            if pvalue_m7m8 >= 0.05:
                g3_pos_sel_var.set("NO")
                g3_pos_sel_entry.config(foreground="black")

        g3_coverage_var.set(coverage.pop(0))
        g3_species_var.set(species.pop(0))

        adapt_less_number = adapt_less.pop(0).split(" ")[0]
        if adapt_less_number == "Not":
            adapt_less_number = 0
        g3_adapt_less_var.set(adapt_less_number)

        adapt_more_number = adapt_more.pop(0).split(" ")[0]
        if adapt_more_number == "Not":
            adapt_more_number = 0
            g3_adapt_more_entry.config(foreground="black")
        elif adapt_more_number == "NA":
            g3_adapt_more_entry.config(foreground="black")
        else:
            g3_adapt_more_entry.config(foreground="magenta")
        g3_adapt_more_var.set(adapt_more_number)

    if gene_name4.get():
        g4_results_var.set(genes.pop(0))
        lrt = lrts.pop(0)
        if lrt == "NA":
            g4_lrt_var.set("NA")
        elif lrt == "(0)":
            g4_lrt_var.set("(0)")
        else:
            g4_lrt_var.set(round(float(lrt), ndigits=4))

        pvalue_m7m8 = pvalues_m7m8.pop(0)
        pvalue_m1m2 = pvalues_m1m2.pop(0)
        if pvalue_m7m8 == "NA":
            g4_pvalue_var.set("NA")
            g4_pos_sel_var.set("NA")
            g4_pos_sel_entry.config(foreground="black")
        else:
            pvalue_m7m8 = float(pvalue_m7m8)
            pvalue_m1m2 = float(pvalue_m1m2)
            if pvalue_m7m8 < 0.05:
                if pvalue_m1m2 < 0.05:
                    g4_pos_sel_var.set("YES")
                else:
                    g4_pos_sel_var.set("(YES)")
                g4_pos_sel_entry.config(foreground="magenta")
            if pvalue_m7m8 <= 0.001:
                g4_pvalue_var.set("<0.001")
            if pvalue_m7m8 > 0.001:
                g4_pvalue_var.set(round(float(pvalue_m7m8), ndigits=4))
            if pvalue_m7m8 >= 0.05:
                g4_pos_sel_var.set("NO")
                g4_pos_sel_entry.config(foreground="black")

        g4_coverage_var.set(coverage.pop(0))
        g4_species_var.set(species.pop(0))

        adapt_less_number = adapt_less.pop(0).split(" ")[0]
        if adapt_less_number == "Not":
            adapt_less_number = 0
        g4_adapt_less_var.set(adapt_less_number)

        adapt_more_number = adapt_more.pop(0).split(" ")[0]
        if adapt_more_number == "Not":
            adapt_more_number = 0
            g4_adapt_more_entry.config(foreground="black")
        elif adapt_more_number == "NA":
            g4_adapt_more_entry.config(foreground="black")
        else:
            g4_adapt_more_entry.config(foreground="magenta")
        g4_adapt_more_var.set(adapt_more_number)

    if gene_name5.get():
        g5_results_var.set(genes.pop(0))
        lrt = lrts.pop(0)
        if lrt == "NA":
            g5_lrt_var.set("NA")
        elif lrt == "(0)":
            g5_lrt_var.set("(0)")
        else:
            g5_lrt_var.set(round(float(lrt), ndigits=4))

        pvalue_m7m8 = pvalues_m7m8.pop(0)
        pvalue_m1m2 = pvalues_m1m2.pop(0)
        if pvalue_m7m8 == "NA":
            g5_pvalue_var.set("NA")
            g5_pos_sel_var.set("NA")
            g5_pos_sel_entry.config(foreground="black")
        else:
            pvalue_m7m8 = float(pvalue_m7m8)
            pvalue_m1m2 = float(pvalue_m1m2)
            if pvalue_m7m8 < 0.05:
                if pvalue_m1m2 < 0.05:
                    g5_pos_sel_var.set("YES")
                else:
                    g5_pos_sel_var.set("(YES)")
                g5_pos_sel_entry.config(foreground="magenta")
            if pvalue_m7m8 <= 0.001:
                g5_pvalue_var.set("<0.001")
            if pvalue_m7m8 > 0.001:
                g5_pvalue_var.set(round(float(pvalue_m7m8), ndigits=4))
            if pvalue_m7m8 >= 0.05:
                g5_pos_sel_var.set("NO")
                g5_pos_sel_entry.config(foreground="black")

        g5_coverage_var.set(coverage.pop(0))
        g5_species_var.set(species.pop(0))

        adapt_less_number = adapt_less.pop(0).split(" ")[0]
        if adapt_less_number == "Not":
            adapt_less_number = 0
        g5_adapt_less_var.set(adapt_less_number)

        adapt_more_number = adapt_more.pop(0).split(" ")[0]
        if adapt_more_number == "Not":
            adapt_more_number = 0
            g5_adapt_more_entry.config(foreground="black")
        elif adapt_more_number == "NA":
            g5_adapt_more_entry.config(foreground="black")
        else:
            g5_adapt_more_entry.config(foreground="magenta")
        g5_adapt_more_var.set(adapt_more_number)


# set up the main window
root = Tk()

root.maxsize(width=1280, height=725)
root.resizable(False, False)

root.title("FREEDA - Finder of Rapidly Evolving Exons in De novo Assemblies")
root.columnconfigure(0, weight=1)
root.columnconfigure(1, weight=1)

# create a mainframe inside the parent frame
mainframe = ttk.Frame(root, padding="2 2 2 2")
mainframe.grid(sticky=(N, W, E, S))
# configure first two columns of the mainframe to behave uniformly
mainframe.columnconfigure(0, weight=1, minsize=550, uniform="group1")
mainframe.columnconfigure(1, weight=1, minsize=700, uniform="group1")
mainframe.rowconfigure(0, weight=1)

# create input frame
input_frame = ttk.Frame(mainframe, relief="sunken", padding="5 5 5 5")
input_frame.grid(column=0, row=0, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
input_frame.columnconfigure(0, weight=3, uniform="group1")
input_frame.columnconfigure((1, 2), weight=1, uniform="group1")
input_frame.columnconfigure(3, weight=2, uniform="group1")

# create logo frame
logo_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
logo_frame.grid(column=0, row=0, columnspan=1, sticky=(N, W, E, S), padx=5, pady=1)

# create settings frame
settings_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
settings_frame.grid(column=1, row=0, columnspan=2, sticky=(N, W, E, S), padx=5, pady=1)
settings_frame.columnconfigure(0, weight=1, uniform="group1")
settings_frame.columnconfigure((1, 2, 3), weight=1, uniform="group1")

# create codon_freq frame
codon_freq_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
codon_freq_frame.grid(column=3, row=0, sticky=(N, W, E, S), padx=5, pady=1)

# create a gene 1 frame
gene1_frame = ttk.Frame(input_frame, relief="ridge", padding="2 2 2 2")
gene1_frame.grid(column=0, row=11, columnspan=4, sticky=(N, W, E, S), padx=5, pady=1)
# let all columns resize
gene1_frame.columnconfigure(0, weight=3, minsize=50, uniform="group1")
gene1_frame.columnconfigure((1, 2), weight=1, minsize=50, uniform="group1")
gene1_frame.columnconfigure(3, weight=2, minsize=50, uniform="group1")

# create a gene 2 frame
gene2_frame = ttk.Frame(input_frame, relief="ridge", padding="2 2 2 2")
gene2_frame.grid(column=0, row=14, columnspan=4, sticky=(N, W, E, S), padx=5, pady=1)
# let all columns resize
gene2_frame.columnconfigure(0, weight=3, uniform="group1")
gene2_frame.columnconfigure((1, 2), weight=1, uniform="group1")
gene2_frame.columnconfigure(3, weight=2, uniform="group1")

# create a gene 3 frame
gene3_frame = ttk.Frame(input_frame, relief="ridge", padding="2 2 2 2")
gene3_frame.grid(column=0, row=17, columnspan=4, sticky=(N, W, E, S), padx=5, pady=1)
# let all columns resize
gene3_frame.columnconfigure(0, weight=3, uniform="group1")
gene3_frame.columnconfigure((1, 2), weight=1, uniform="group1")
gene3_frame.columnconfigure(3, weight=2, uniform="group1")

# create a gene 4 frame
gene4_frame = ttk.Frame(input_frame, relief="ridge", padding="2 2 2 2")
gene4_frame.grid(column=0, row=18, columnspan=4, sticky=(N, W, E, S), padx=5, pady=1)
# let all columns resize
gene4_frame.columnconfigure(0, weight=3, uniform="group1")
gene4_frame.columnconfigure((1, 2), weight=1, uniform="group1")
gene4_frame.columnconfigure(3, weight=2, uniform="group1")

# create a gene 5 frame
gene5_frame = ttk.Frame(input_frame, relief="ridge", padding="2 2 2 2")
gene5_frame.grid(column=0, row=19, columnspan=4, sticky=(N, W, E, S), padx=5, pady=1)
# let all columns resize
gene5_frame.columnconfigure(0, weight=3, uniform="group1")
gene5_frame.columnconfigure((1, 2), weight=1, uniform="group1")
gene5_frame.columnconfigure(3, weight=2, uniform="group1")

# create a gene working directory frame
wdir_frame = ttk.Frame(input_frame, relief="ridge", padding="2 2 2 2")
wdir_frame.grid(column=0, row=20, columnspan=5, sticky=(N, W, E, S), padx=5, pady=1)
# let all columns resize
wdir_frame.columnconfigure(0, weight=3, uniform="group1")
wdir_frame.columnconfigure(1, weight=5, uniform="group1")
wdir_frame.columnconfigure((2, 3), weight=2, uniform="group1")

# create exclude frame
exclude_frame = ttk.Frame(input_frame, relief="ridge", padding="2 2 2 2")
exclude_frame.grid(column=0, row=21, columnspan=5, sticky=(N, W, E, S), padx=5, pady=1)
# let all columns resize
exclude_frame.columnconfigure(0, weight=3, uniform="group1")
exclude_frame.columnconfigure(1, weight=5, uniform="group1")
exclude_frame.columnconfigure((2, 3), weight=2, uniform="group1")

# create output frame
output_frame = ttk.Frame(mainframe, relief="sunken", padding="5 5 5 5")
output_frame.grid(column=1, row=0, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
output_frame.columnconfigure((0, 1, 2, 3, 4, 5), weight=1, uniform="group1")

# create a logging window
logging_frame = ttk.Frame(output_frame, relief="ridge", padding="5 5 5 5")
logging_frame.grid(column=0, row=0, columnspan=9, sticky=(N, W, E, S), padx=5, pady=5)

# create results frame
results_labelframe = ttk.LabelFrame(output_frame, text="Results window (output only for F3X4)")
results_labelframe.grid(column=0, row=1, columnspan=8, padx=5, pady=2, sticky=(N, W))
results_labelframe.columnconfigure((0, 1, 2, 3, 4, 5, 6, 7), weight=1, uniform="group1")

# create results labels frame upper
results_labels_upper = ttk.Frame(results_labelframe, padding="2 2 2 2")
results_labels_upper.grid(column=0, row=2, columnspan=8, sticky=(N, W, E, S), padx=5, pady=2)
results_labels_upper.columnconfigure((0, 1, 2, 3, 4, 5, 6, 7), weight=1, uniform="group1")

# create results labels frame
results_labels = ttk.Frame(results_labelframe, padding="2 2 2 2")
results_labels.grid(column=0, row=3, columnspan=8, sticky=(N, W, E, S), padx=5, pady=2)
results_labels.columnconfigure((0, 1, 2, 3, 4, 5, 6, 7), weight=1, uniform="group1")

# create user gene 1 result frame
g1_results_frame = ttk.Frame(results_labelframe, relief="ridge", padding="5 5 5 5")
g1_results_frame.grid(column=0, row=4, columnspan=8, sticky=(N, W, E, S), padx=5, pady=5)
g1_results_frame.columnconfigure((0, 1, 2, 3, 4, 5, 6, 7), weight=1, uniform="group1")

# create user gene 2 result frame
g2_results_frame = ttk.Frame(results_labelframe, relief="ridge", padding="5 5 5 5")
g2_results_frame.grid(column=0, row=5, columnspan=8, sticky=(N, W, E, S), padx=5, pady=5)
g2_results_frame.columnconfigure((0, 1, 2, 3, 4, 5, 6, 7), weight=1, uniform="group1")

# create user gene 3 result frame
g3_results_frame = ttk.Frame(results_labelframe, relief="ridge", padding="5 5 5 5")
g3_results_frame.grid(column=0, row=6, columnspan=8, sticky=(N, W, E, S), padx=5, pady=5)
g3_results_frame.columnconfigure((0, 1, 2, 3, 4, 5, 6, 7), weight=1, uniform="group1")

# create user gene 4 result frame
g4_results_frame = ttk.Frame(results_labelframe, relief="ridge", padding="5 5 5 5")
g4_results_frame.grid(column=0, row=7, columnspan=8, sticky=(N, W, E, S), padx=5, pady=5)
g4_results_frame.columnconfigure((0, 1, 2, 3, 4, 5, 6, 7), weight=1, uniform="group1")

# create user gene 5 result frame
g5_results_frame = ttk.Frame(results_labelframe, relief="ridge", padding="5 5 5 5")
g5_results_frame.grid(column=0, row=8, columnspan=8, sticky=(N, W, E, S), padx=5, pady=5)
g5_results_frame.columnconfigure((0, 1, 2, 3, 4, 5, 6, 7), weight=1, uniform="group1")

# CHECKS
check_gene_name_wrapper = (settings_frame.register(check_gene_name), "%P", "%V")
check_functional_residues = (settings_frame.register(check_functional_residues), "%P", "%V")
check_label = (settings_frame.register(check_label), "%P", "%V")
check_excluded_species = (settings_frame.register(check_excluded_species), "%P", "%V")

# ERRORS
error_message1 = StringVar()
message1 = "Invalid gene name (follow pattern: Cenpo for mouse and : CENPO for others)"
error_message2 = StringVar()
message2 = "Invalid residue number (follow pattern: 100)"
# error labels
error_label1 = ttk.Label(input_frame, font="TkSmallCaptionFont", foreground="magenta", textvariable=error_message1)
error_label1.grid(column=0, row=22, columnspan=4, padx=5, pady=1, sticky=(W))
error_label2 = ttk.Label(input_frame, font="TkSmallCaptionFont", foreground="magenta", textvariable=error_message2)
error_label2.grid(column=0, row=23, columnspan=4, padx=5, pady=1, sticky=(W))

# LOGO
canvas = Canvas(logo_frame, width=200, height=50)
canvas.pack(expand=True)
freeda_img = ImageTk.PhotoImage(Image.open(pyinstaller_compatibility.resource_path("freeda_logo.png")))
canvas.create_image(100, 30, anchor=CENTER, image=freeda_img)

# LOGGER
raise_logger()

# CLADE
clade = StringVar()
#ttk.Label(settings_frame, text="Clade").grid(column=0, row=0, pady=5, sticky=(W))
rodents = ttk.Radiobutton(settings_frame, text="Mouse (Murinae)", variable=clade, value="Mm")
rodents.grid(column=0, row=1, columnspan=3, sticky=(W))
primates = ttk.Radiobutton(settings_frame, text="Human (Primates)", variable=clade, value="Hs")
primates.grid(column=0, row=2, columnspan=4, sticky=(W))
carnivores = ttk.Radiobutton(settings_frame, text="Dog (Carnivora)", variable=clade, value="Cf")
carnivores.grid(column=0, row=3, columnspan=3, sticky=(W))
birds = ttk.Radiobutton(settings_frame, text="Chicken (Phasian.)", variable=clade, value="Gg")
birds.grid(column=0, row=4, columnspan=3, sticky=(W))

# CODON FREQUENCY
ttk.Label(codon_freq_frame, text="Codon frequencies").grid(column=0, row=0, sticky=(W))
codon_freq = StringVar()
# set F3X4 codon frequency as default
codon_freq.set("F3X4")
codon_freq_F3X4 = ttk.Radiobutton(codon_freq_frame, text="F3X4 (default)", variable=codon_freq, value="F3X4")
codon_freq_F3X4.grid(column=0, row=1, sticky=(W))
codon_freq_F61 = ttk.Radiobutton(codon_freq_frame, text="F3X4 and F61", variable=codon_freq, value="F3X4,F61")
codon_freq_F61.grid(column=0, row=2, sticky=(W))

# RESIDUES
ttk.Label(input_frame, text="Indicate functional residues -----> ").grid(column=0, row=9, columnspan=2,
                                                                     sticky=(W, E), padx=5)
ttk.Label(input_frame, text="start").grid(column=1, row=9)
ttk.Label(input_frame, text="end").grid(column=2, row=9)
ttk.Label(input_frame, text="label").grid(column=3, row=9)

# USER INPUT GENE 1
gene_name1 = StringVar()
ttk.Label(gene1_frame, text="Gene name").grid(column=0, row=11, padx=6, pady=1, sticky=(W))
name1 = ttk.Entry(gene1_frame, textvariable=gene_name1, validate="all", validatecommand=check_gene_name_wrapper)
name1.grid(column=0, row=12, padx=5, pady=1, sticky=(W))
dup1_var = BooleanVar()
dup1_button = ttk.Checkbutton(gene1_frame, text="Duplication expected", variable=dup1_var, onvalue=1, offvalue=0)
dup1_button.grid(column=0, row=13, padx=6, pady=1, sticky=(W))

s11_start = StringVar()
s11_end = StringVar()
s11_label = StringVar()
site11_start = ttk.Entry(gene1_frame, textvariable=s11_start, validate="all", validatecommand=check_functional_residues)
site11_start.grid(column=1, row=11, padx=5, pady=1, sticky=(W))
site11_end = ttk.Entry(gene1_frame, textvariable=s11_end, validate="all", validatecommand=check_functional_residues)
site11_end.grid(column=2, row=11, padx=5, pady=1, sticky=(W))
site11_label = ttk.Entry(gene1_frame, textvariable=s11_label, validate="all", validatecommand=check_label)
site11_label.grid(column=3, row=11, padx=5, pady=1, sticky=(W))

s12_start = StringVar()
s12_end = StringVar()
s12_label = StringVar()
site12_start = ttk.Entry(gene1_frame, textvariable=s12_start, validate="all", validatecommand=check_functional_residues)
site12_start.grid(column=1, row=12, padx=5, pady=1, sticky=(W))
site12_end = ttk.Entry(gene1_frame, textvariable=s12_end, validate="all", validatecommand=check_functional_residues)
site12_end.grid(column=2, row=12, padx=5, pady=1, sticky=(W))
site12_label = ttk.Entry(gene1_frame, textvariable=s12_label, validate="all", validatecommand=check_label)
site12_label.grid(column=3, row=12, padx=5, pady=1, sticky=(W))

s13_start = StringVar()
s13_end = StringVar()
s13_label = StringVar()
site13_start = ttk.Entry(gene1_frame, textvariable=s13_start, validate="all", validatecommand=check_functional_residues)
site13_start.grid(column=1, row=13, padx=5, pady=1, sticky=(W))
site13_end = ttk.Entry(gene1_frame, textvariable=s13_end, validate="all", validatecommand=check_functional_residues)
site13_end.grid(column=2, row=13, padx=5, pady=1, sticky=(W))
site13_label = ttk.Entry(gene1_frame, textvariable=s13_label, validate="all", validatecommand=check_label)
site13_label.grid(column=3, row=13, padx=5, pady=1, sticky=(W))

# USER INPUT GENE 2
gene_name2 = StringVar()
ttk.Label(gene2_frame, text="Gene name").grid(column=0, row=14, padx=6, pady=1, sticky=(W))
name2 = ttk.Entry(gene2_frame, textvariable=gene_name2, validate="all", validatecommand=check_gene_name_wrapper)
name2.grid(column=0, row=15, padx=5, pady=1, sticky=(W))
dup2_var = BooleanVar()
dup2_button = ttk.Checkbutton(gene2_frame, text="Duplication expected", variable=dup2_var, onvalue=1, offvalue=0)
dup2_button.grid(column=0, row=16, padx=6, pady=1, sticky=(W))

s21_start = StringVar()
s21_end = StringVar()
s21_label = StringVar()
site21_start = ttk.Entry(gene2_frame, textvariable=s21_start, validate="all", validatecommand=check_functional_residues)
site21_start.grid(column=1, row=14, padx=5, pady=1, sticky=(W))
site21_end = ttk.Entry(gene2_frame, textvariable=s21_end, validate="all", validatecommand=check_functional_residues)
site21_end.grid(column=2, row=14, padx=5, pady=1, sticky=(W))
site21_label = ttk.Entry(gene2_frame, textvariable=s21_label, validate="all", validatecommand=check_label)
site21_label.grid(column=3, row=14, padx=5, pady=1, sticky=(W))

s22_start = StringVar()
s22_end = StringVar()
s22_label = StringVar()
site22_start = ttk.Entry(gene2_frame, textvariable=s22_start, validate="all", validatecommand=check_functional_residues)
site22_start.grid(column=1, row=15, padx=5, pady=1, sticky=(W))
site22_end = ttk.Entry(gene2_frame, textvariable=s22_end, validate="all", validatecommand=check_functional_residues)
site22_end.grid(column=2, row=15, padx=5, pady=1, sticky=(W))
site22_label = ttk.Entry(gene2_frame, textvariable=s22_label, validate="all", validatecommand=check_label)
site22_label.grid(column=3, row=15, padx=5, pady=1, sticky=(W))

s23_start = StringVar()
s23_end = StringVar()
s23_label = StringVar()
site23_start = ttk.Entry(gene2_frame, textvariable=s23_start, validate="all", validatecommand=check_functional_residues)
site23_start.grid(column=1, row=16, padx=5, pady=1, sticky=(W))
site23_end = ttk.Entry(gene2_frame, textvariable=s23_end, validate="all", validatecommand=check_functional_residues)
site23_end.grid(column=2, row=16, padx=5, pady=1, sticky=(W))
site23_label = ttk.Entry(gene2_frame, textvariable=s23_label, validate="all", validatecommand=check_label)
site23_label.grid(column=3, row=16, padx=5, pady=1, sticky=(W))

# USER INPUT GENE 3
gene_name3 = StringVar()
ttk.Label(gene3_frame, text="Gene name").grid(column=0, row=17, padx=6, pady=1, sticky=(W))
name3 = ttk.Entry(gene3_frame, textvariable=gene_name3, validate="all", validatecommand=check_gene_name_wrapper)
name3.grid(column=0, row=18, padx=5, pady=1, sticky=(W))
dup3_var = BooleanVar()
dup3_button = ttk.Checkbutton(gene3_frame, text="Duplication expected", variable=dup3_var, onvalue=1, offvalue=0)
dup3_button.grid(column=0, row=19, padx=6, pady=1, sticky=(W))

s31_start = StringVar()
s31_end = StringVar()
s31_label = StringVar()
site31_start = ttk.Entry(gene3_frame, textvariable=s31_start, validate="all", validatecommand=check_functional_residues)
site31_start.grid(column=1, row=17, padx=5, pady=1, sticky=(W))
site31_end = ttk.Entry(gene3_frame, textvariable=s31_end, validate="all", validatecommand=check_functional_residues)
site31_end.grid(column=2, row=17, padx=5, pady=1, sticky=(W))
site31_label = ttk.Entry(gene3_frame, textvariable=s31_label, validate="all", validatecommand=check_label)
site31_label.grid(column=3, row=17, padx=5, pady=1, sticky=(W))

s32_start = StringVar()
s32_end = StringVar()
s32_label = StringVar()
site32_start = ttk.Entry(gene3_frame, textvariable=s32_start, validate="all", validatecommand=check_functional_residues)
site32_start.grid(column=1, row=18, padx=5, pady=1, sticky=(W))
site32_end = ttk.Entry(gene3_frame, textvariable=s32_end, validate="all", validatecommand=check_functional_residues)
site32_end.grid(column=2, row=18, padx=5, pady=1, sticky=(W))
site32_label = ttk.Entry(gene3_frame, textvariable=s32_label, validate="all", validatecommand=check_label)
site32_label.grid(column=3, row=18, padx=5, pady=1, sticky=(W))

s33_start = StringVar()
s33_end = StringVar()
s33_label = StringVar()
site33_start = ttk.Entry(gene3_frame, textvariable=s33_start, validate="all", validatecommand=check_functional_residues)
site33_start.grid(column=1, row=19, padx=5, pady=1, sticky=(W))
site33_end = ttk.Entry(gene3_frame, textvariable=s33_end, validate="all", validatecommand=check_functional_residues)
site33_end.grid(column=2, row=19, padx=5, pady=1, sticky=(W))
site33_label = ttk.Entry(gene3_frame, textvariable=s33_label, validate="all", validatecommand=check_label)
site33_label.grid(column=3, row=19, padx=5, pady=1, sticky=(W))

# USER INPUT GENE 4
gene_name4 = StringVar()
ttk.Label(gene4_frame, text="Gene name").grid(column=0, row=20, padx=6, pady=1, sticky=(W))
name4 = ttk.Entry(gene4_frame, textvariable=gene_name4, validate="all", validatecommand=check_gene_name_wrapper)
name4.grid(column=0, row=21, padx=5, pady=1, sticky=(W))
dup4_var = BooleanVar()
dup4_button = ttk.Checkbutton(gene4_frame, text="Duplication expected", variable=dup4_var, onvalue=1, offvalue=0)
dup4_button.grid(column=0, row=22, padx=6, pady=1, sticky=(W))

s41_start = StringVar()
s41_end = StringVar()
s41_label = StringVar()
site41_start = ttk.Entry(gene4_frame, textvariable=s41_start, validate="all", validatecommand=check_functional_residues)
site41_start.grid(column=1, row=20, padx=5, pady=1, sticky=(W))
site41_end = ttk.Entry(gene4_frame, textvariable=s41_end, validate="all", validatecommand=check_functional_residues)
site41_end.grid(column=2, row=20, padx=5, pady=1, sticky=(W))
site41_label = ttk.Entry(gene4_frame, textvariable=s41_label, validate="all", validatecommand=check_label)
site41_label.grid(column=3, row=20, padx=5, pady=1, sticky=(W))

s42_start = StringVar()
s42_end = StringVar()
s42_label = StringVar()
site42_start = ttk.Entry(gene4_frame, textvariable=s42_start, validate="all", validatecommand=check_functional_residues)
site42_start.grid(column=1, row=21, padx=5, pady=1, sticky=(W))
site42_end = ttk.Entry(gene4_frame, textvariable=s42_end, validate="all", validatecommand=check_functional_residues)
site42_end.grid(column=2, row=21, padx=5, pady=1, sticky=(W))
site42_label = ttk.Entry(gene4_frame, textvariable=s42_label, validate="all", validatecommand=check_label)
site42_label.grid(column=3, row=21, padx=5, pady=1, sticky=(W))

s43_start = StringVar()
s43_end = StringVar()
s43_label = StringVar()
site43_start = ttk.Entry(gene4_frame, textvariable=s43_start, validate="all", validatecommand=check_functional_residues)
site43_start.grid(column=1, row=22, padx=5, pady=1, sticky=(W))
site43_end = ttk.Entry(gene4_frame, textvariable=s43_end, validate="all", validatecommand=check_functional_residues)
site43_end.grid(column=2, row=22, padx=5, pady=1, sticky=(W))
site43_label = ttk.Entry(gene4_frame, textvariable=s43_label, validate="all", validatecommand=check_label)
site43_label.grid(column=3, row=22, padx=5, pady=1, sticky=(W))

# USER INPUT GENE 5
gene_name5 = StringVar()
ttk.Label(gene5_frame, text="Gene name").grid(column=0, row=23, padx=6, pady=1, sticky=(W))
name5 = ttk.Entry(gene5_frame, textvariable=gene_name5, validate="all", validatecommand=check_gene_name_wrapper)
name5.grid(column=0, row=24, padx=5, pady=1, sticky=(W))
dup5_var = BooleanVar()
dup5_button = ttk.Checkbutton(gene5_frame, text="Duplication expected", variable=dup5_var, onvalue=1, offvalue=0)
dup5_button.grid(column=0, row=25, padx=6, pady=1, sticky=(W))

s51_start = StringVar()
s51_end = StringVar()
s51_label = StringVar()
site51_start = ttk.Entry(gene5_frame, textvariable=s51_start, validate="all", validatecommand=check_functional_residues)
site51_start.grid(column=1, row=23, padx=5, pady=1, sticky=(W))
site51_end = ttk.Entry(gene5_frame, textvariable=s51_end, validate="all", validatecommand=check_functional_residues)
site51_end.grid(column=2, row=23, padx=5, pady=1, sticky=(W))
site51_label = ttk.Entry(gene5_frame, textvariable=s51_label, validate="all", validatecommand=check_label)
site51_label.grid(column=3, row=23, padx=5, pady=1, sticky=(W))

s52_start = StringVar()
s52_end = StringVar()
s52_label = StringVar()
site52_start = ttk.Entry(gene5_frame, textvariable=s52_start, validate="all", validatecommand=check_functional_residues)
site52_start.grid(column=1, row=24, padx=5, pady=1, sticky=(W))
site52_end = ttk.Entry(gene5_frame, textvariable=s52_end, validate="all", validatecommand=check_functional_residues)
site52_end.grid(column=2, row=24, padx=5, pady=1, sticky=(W))
site52_label = ttk.Entry(gene5_frame, textvariable=s52_label, validate="all", validatecommand=check_label)
site52_label.grid(column=3, row=24, padx=5, pady=1, sticky=(W))

s53_start = StringVar()
s53_end = StringVar()
s53_label = StringVar()
site53_start = ttk.Entry(gene5_frame, textvariable=s53_start, validate="all", validatecommand=check_functional_residues)
site53_start.grid(column=1, row=25, padx=5, pady=1, sticky=(W))
site53_end = ttk.Entry(gene5_frame, textvariable=s53_end, validate="all", validatecommand=check_functional_residues)
site53_end.grid(column=2, row=25, padx=5, pady=1, sticky=(W))
site53_label = ttk.Entry(gene5_frame, textvariable=s53_label, validate="all", validatecommand=check_label)
site53_label.grid(column=3, row=25, padx=5, pady=1, sticky=(W))

# EXCLUDE SPECIES
exclude_species_var = StringVar()
exclude_species_entry = ttk.Entry(exclude_frame, textvariable=exclude_species_var, validate="all",
                                  validatecommand=check_excluded_species)
exclude_species_entry.grid(column=1, row=21, padx=5, pady=1, sticky=(W, E))
ttk.Label(exclude_frame, text="   Exclude species: ").grid(column=0, row=21, padx=5, sticky=(W))
ttk.Label(exclude_frame, text="(optional; e.g. Ha Gs)").grid(column=2, row=21, columnspan=2, sticky=(W))

# BUTTONS
result_path_var = StringVar()

# Working directory button
wdir_button = ttk.Button(wdir_frame, text="Set directory", command=get_wdir)
wdir_button.grid(column=0, row=0, sticky=(N, W, E, S), padx=1, pady=1)
wdirectory = StringVar()
wdir_entry = ttk.Entry(wdir_frame, textvariable=wdirectory)
wdir_entry.grid(column=1, row=0, sticky=(N, W, E, S), padx=5, pady=1)  # columnspan=4

# ANALYZE button
analyze_button = ttk.Button(wdir_frame, text="Analyze", state="normal", command=thread_freeda)  # default="active"
analyze_button.grid(column=2, row=0, pady=1, sticky=(W))  # columnspan=2

# ABORT button
abort_button = ttk.Button(wdir_frame, text="ABORT", state="normal", command=abort_freeda)  # default="active"
abort_button.grid(column=3, row=0,  pady=1, sticky=(W))

# LABELS
ttk.Label(results_labels_upper, text="nr of adapt. residues").grid(column=7, row=2, columnspan=2, sticky=(E), padx=15)
ttk.Label(results_labels, text="Gene").grid(column=0, row=3)
ttk.Label(results_labels, text="Pos. select.").grid(column=1, row=3)
ttk.Label(results_labels, text="LRT").grid(column=2, row=3)
ttk.Label(results_labels, text="p-value").grid(column=3, row=3)
ttk.Label(results_labels, text="CDS cover.").grid(column=4, row=3)
ttk.Label(results_labels, text="species").grid(column=5, row=3)
ttk.Label(results_labels, text="pr < 0.9").grid(column=6, row=3)
ttk.Label(results_labels, text="pr >= 0.9").grid(column=7, row=3)

# ENTRIES
# gene 1 results
g1_results_var = StringVar()
g1_results_entry = ttk.Entry(g1_results_frame, state="disabled", text=g1_results_var, justify='center')
g1_results_entry.config(foreground="black")  # text will be black despite disabled state
g1_results_entry.grid(column=0, row=0, sticky=(W))
g1_pos_sel_var = StringVar()
g1_pos_sel_entry = ttk.Entry(g1_results_frame, state="disabled", text=g1_pos_sel_var, justify='center')
g1_pos_sel_entry.grid(column=1, row=0, sticky=(W))
g1_pos_sel_entry.config(foreground="black")  # text will be black despite disabled state
g1_lrt_var = StringVar()
g1_lrt_entry = ttk.Entry(g1_results_frame, state="disabled", text=g1_lrt_var, justify='center')
g1_lrt_entry.grid(column=2, row=0, sticky=(W))
g1_lrt_entry.config(foreground="black")  # text will be black despite disabled state
g1_pvalue_var = StringVar()
g1_pvalue_entry = ttk.Entry(g1_results_frame, state="disabled", text=g1_pvalue_var, justify='center')
g1_pvalue_entry.grid(column=3, row=0, sticky=(W))
g1_pvalue_entry.config(foreground="black")  # text will be black despite disabled state
g1_coverage_var = StringVar()
g1_coverage_entry = ttk.Entry(g1_results_frame, state="disabled", text=g1_coverage_var, justify='center')
g1_coverage_entry.grid(column=4, row=0, sticky=(W))
g1_coverage_entry.config(foreground="black")  # text will be black despite disabled state
g1_species_var = StringVar()
g1_species_entry = ttk.Entry(g1_results_frame, state="disabled", text=g1_species_var, justify='center')
g1_species_entry.grid(column=5, row=0, sticky=(W))
g1_species_entry.config(foreground="black")  # text will be black despite disabled state
g1_adapt_less_var = StringVar()
g1_adapt_less_entry = ttk.Entry(g1_results_frame, state="disabled", text=g1_adapt_less_var, justify='center')
g1_adapt_less_entry.grid(column=6, row=0, sticky=(W))
g1_adapt_less_entry.config(foreground="black")  # text will be black despite disabled state
g1_adapt_more_var = StringVar()
g1_adapt_more_entry = ttk.Entry(g1_results_frame, state="disabled", text=g1_adapt_more_var, justify='center')
g1_adapt_more_entry.grid(column=7, row=0, sticky=(W))
g1_adapt_more_entry.config(foreground="black")  # text will be black despite disabled state


# gene 2 results
g2_results_var = StringVar()
g2_results_entry = ttk.Entry(g2_results_frame, state="disabled", text=g2_results_var, justify='center')
g2_results_entry.grid(column=0, row=0, sticky=(W))
g2_results_entry.config(foreground="black")  # text will be black despite disabled state
g2_pos_sel_var = StringVar()
g2_pos_sel_entry = ttk.Entry(g2_results_frame, state="disabled", text=g2_pos_sel_var, justify='center')
g2_pos_sel_entry.grid(column=1, row=0, sticky=(W))
g2_pos_sel_entry.config(foreground="black")  # text will be black despite disabled state
g2_lrt_var = StringVar()
g2_lrt_entry = ttk.Entry(g2_results_frame, state="disabled", text=g2_lrt_var, justify='center')
g2_lrt_entry.grid(column=2, row=0, sticky=(W))
g2_lrt_entry.config(foreground="black")  # text will be black despite disabled state
g2_pvalue_var = StringVar()
g2_pvalue_entry = ttk.Entry(g2_results_frame, state="disabled", text=g2_pvalue_var, justify='center')
g2_pvalue_entry.grid(column=3, row=0, sticky=(W))
g2_pvalue_entry.config(foreground="black")  # text will be black despite disabled state
g2_coverage_var = StringVar()
g2_coverage_entry = ttk.Entry(g2_results_frame, state="disabled", text=g2_coverage_var, justify='center')
g2_coverage_entry.grid(column=4, row=0, sticky=(W))
g2_coverage_entry.config(foreground="black")  # text will be black despite disabled state
g2_species_var = StringVar()
g2_species_entry = ttk.Entry(g2_results_frame, state="disabled", text=g2_species_var, justify='center')
g2_species_entry.grid(column=5, row=0, sticky=(W))
g2_species_entry.config(foreground="black")  # text will be black despite disabled state
g2_adapt_less_var = StringVar()
g2_adapt_less_entry = ttk.Entry(g2_results_frame, state="disabled", text=g2_adapt_less_var, justify='center')
g2_adapt_less_entry.grid(column=6, row=0, sticky=(W))
g2_adapt_less_entry.config(foreground="black")  # text will be black despite disabled state
g2_adapt_more_var = StringVar()
g2_adapt_more_entry = ttk.Entry(g2_results_frame, state="disabled", text=g2_adapt_more_var, justify='center')
g2_adapt_more_entry.grid(column=7, row=0, sticky=(W))
g2_adapt_more_entry.config(foreground="black")  # text will be black despite disabled state

# gene 3 results
g3_results_var = StringVar()
g3_results_entry = ttk.Entry(g3_results_frame, state="disabled", text=g3_results_var, justify='center')
g3_results_entry.grid(column=0, row=0, sticky=(W))
g3_results_entry.config(foreground="black")  # text will be black despite disabled state
g3_pos_sel_var = StringVar()
g3_pos_sel_entry = ttk.Entry(g3_results_frame, state="disabled", text=g3_pos_sel_var, justify='center')
g3_pos_sel_entry.grid(column=1, row=0, sticky=(W))
g3_pos_sel_entry.config(foreground="black")  # text will be black despite disabled state
g3_lrt_var = StringVar()
g3_lrt_entry = ttk.Entry(g3_results_frame, state="disabled", text=g3_lrt_var, justify='center')
g3_lrt_entry.grid(column=2, row=0, sticky=(W))
g3_lrt_entry.config(foreground="black")  # text will be black despite disabled state
g3_pvalue_var = StringVar()
g3_pvalue_entry = ttk.Entry(g3_results_frame, state="disabled", text=g3_pvalue_var, justify='center')
g3_pvalue_entry.grid(column=3, row=0, sticky=(W))
g3_pvalue_entry.config(foreground="black")  # text will be black despite disabled state
g3_coverage_var = StringVar()
g3_coverage_entry = ttk.Entry(g3_results_frame, state="disabled", text=g3_coverage_var, justify='center')
g3_coverage_entry.grid(column=4, row=0, sticky=(W))
g3_coverage_entry.config(foreground="black")  # text will be black despite disabled state
g3_species_var = StringVar()
g3_species_entry = ttk.Entry(g3_results_frame, state="disabled", text=g3_species_var, justify='center')
g3_species_entry.grid(column=5, row=0, sticky=(W))
g3_species_entry.config(foreground="black")  # text will be black despite disabled state
g3_adapt_less_var = StringVar()
g3_adapt_less_entry = ttk.Entry(g3_results_frame, state="disabled", text=g3_adapt_less_var, justify='center')
g3_adapt_less_entry.grid(column=6, row=0, sticky=(W))
g3_adapt_less_entry.config(foreground="black")  # text will be black despite disabled state
g3_adapt_more_var = StringVar()
g3_adapt_more_entry = ttk.Entry(g3_results_frame, state="disabled", text=g3_adapt_more_var, justify='center')
g3_adapt_more_entry.grid(column=7, row=0, sticky=(W))
g3_adapt_more_entry.config(foreground="black")  # text will be black despite disabled state


# gene 4 results
g4_results_var = StringVar()
g4_results_entry = ttk.Entry(g4_results_frame, state="disabled", text=g4_results_var, justify='center')
g4_results_entry.grid(column=0, row=0, sticky=(W))
g4_results_entry.config(foreground="black")  # text will be black despite disabled state
g4_pos_sel_var = StringVar()
g4_pos_sel_entry = ttk.Entry(g4_results_frame, state="disabled", text=g4_pos_sel_var, justify='center')
g4_pos_sel_entry.grid(column=1, row=0, sticky=(W))
g4_pos_sel_entry.config(foreground="black")  # text will be black despite disabled state
g4_lrt_var = StringVar()
g4_lrt_entry = ttk.Entry(g4_results_frame, state="disabled", text=g4_lrt_var, justify='center')
g4_lrt_entry.grid(column=2, row=0, sticky=(W))
g4_lrt_entry.config(foreground="black")  # text will be black despite disabled state
g4_pvalue_var = StringVar()
g4_pvalue_entry = ttk.Entry(g4_results_frame, state="disabled", text=g4_pvalue_var, justify='center')
g4_pvalue_entry.grid(column=3, row=0, sticky=(W))
g4_pvalue_entry.config(foreground="black")  # text will be black despite disabled state
g4_coverage_var = StringVar()
g4_coverage_entry = ttk.Entry(g4_results_frame, state="disabled", text=g4_coverage_var, justify='center')
g4_coverage_entry.grid(column=4, row=0, sticky=(W))
g4_coverage_entry.config(foreground="black")  # text will be black despite disabled state
g4_species_var = StringVar()
g4_species_entry = ttk.Entry(g4_results_frame, state="disabled", text=g4_species_var, justify='center')
g4_species_entry.grid(column=5, row=0, sticky=(W))
g4_species_entry.config(foreground="black")  # text will be black despite disabled state
g4_adapt_less_var = StringVar()
g4_adapt_less_entry = ttk.Entry(g4_results_frame, state="disabled", text=g4_adapt_less_var, justify='center')
g4_adapt_less_entry.grid(column=6, row=0, sticky=(W))
g4_adapt_less_entry.config(foreground="black")  # text will be black despite disabled state
g4_adapt_more_var = StringVar()
g4_adapt_more_entry = ttk.Entry(g4_results_frame, state="disabled", text=g4_adapt_more_var, justify='center')
g4_adapt_more_entry.grid(column=7, row=0, sticky=(W))
g4_adapt_more_entry.config(foreground="black")  # text will be black despite disabled state


# gene 5 results
g5_results_var = StringVar()
g5_results_entry = ttk.Entry(g5_results_frame, state="disabled", text=g5_results_var, justify='center')
g5_results_entry.grid(column=0, row=0, sticky=(W))
g5_results_entry.config(foreground="black")  # text will be black despite disabled state
g5_pos_sel_var = StringVar()
g5_pos_sel_entry = ttk.Entry(g5_results_frame, state="disabled", text=g5_pos_sel_var, justify='center')
g5_pos_sel_entry.grid(column=1, row=0, sticky=(W))
g5_pos_sel_entry.config(foreground="black")  # text will be black despite disabled state
g5_lrt_var = StringVar()
g5_lrt_entry = ttk.Entry(g5_results_frame, state="disabled", text=g5_lrt_var, justify='center')
g5_lrt_entry.grid(column=2, row=0, sticky=(W))
g5_lrt_entry.config(foreground="black")  # text will be black despite disabled state
g5_pvalue_var = StringVar()
g5_pvalue_entry = ttk.Entry(g5_results_frame, state="disabled", text=g5_pvalue_var, justify='center')
g5_pvalue_entry.grid(column=3, row=0, sticky=(W))
g5_pvalue_entry.config(foreground="black")  # text will be black despite disabled state
g5_coverage_var = StringVar()
g5_coverage_entry = ttk.Entry(g5_results_frame, state="disabled", text=g5_coverage_var, justify='center')
g5_coverage_entry.grid(column=4, row=0, sticky=(W))
g5_coverage_entry.config(foreground="black")  # text will be black despite disabled state
g5_species_var = StringVar()
g5_species_entry = ttk.Entry(g5_results_frame, state="disabled", text=g5_species_var, justify='center')
g5_species_entry.grid(column=5, row=0, sticky=(W))
g5_species_entry.config(foreground="black")  # text will be black despite disabled state
g5_adapt_less_var = StringVar()
g5_adapt_less_entry = ttk.Entry(g5_results_frame, state="disabled", text=g5_adapt_less_var, justify='center')
g5_adapt_less_entry.grid(column=6, row=0, sticky=(W))
g5_adapt_less_entry.config(foreground="black")  # text will be black despite disabled state
g5_adapt_more_var = StringVar()
g5_adapt_more_entry = ttk.Entry(g5_results_frame, state="disabled", text=g5_adapt_more_var, justify='center')
g5_adapt_more_entry.grid(column=7, row=0, sticky=(W))
g5_adapt_more_entry.config(foreground="black")  # text will be black despite disabled state


try:
    root.mainloop()
except KeyboardInterrupt as kie:
    print("exception in main", kie)
