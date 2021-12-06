#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:36:16 2021

@author: damian
Main module of the freeda package. Takes user input and performs automatic input extraction, tblastn, exon finding
and molecular evolution analysis (PAML) followed by overlay of putative adaptive sites onto 3D structure (PyMOL).

"""

# TODO
#       0) Finish logging -> DONE
#       1) "/Users/damian/PycharmProjects/freeda_2.0/freeda/paml_visualizer.py:844:
#       UserWarning: Starting a Matplotlib GUI outside of the main thread will likely fail. -> FIXED by "use("agg")
#       2) STOP button -> DONE
#       3) Analyse button should block after pushing -> DONE
#       4) Allow deleting the first letter -> DONE
#       5) Connect sites and duplicated expected
#       6) Make result window -> DONE
#       7) Make the main window smaller
#       8) Make Analyse button non-responsive if no wdir, clade, threshold or gene name -> DONE


from freeda import input_extractor
from freeda import tblastn
from freeda import exon_extractor
from freeda import paml_launcher
from freeda import paml_visualizer
from freeda import structure_builder
from freeda import genomes_preprocessing
from freeda import TextHandler
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from PIL import ImageTk, Image
import tkinter.scrolledtext as ScrolledText
import os
import re
import logging
import threading


def thread_freeda():
    """Runs freeda main function in another thread"""

    if check_input():
        freeda_thread = threading.Thread(target=freeda_pipeline)
        # makes a deamon thread allowing it to be killed anytime (though not advised generally)
        freeda_thread.daemon = True
        freeda_thread.start()
    else:
        logging.info("\n--------------- TRY AGAIN :D -------------------\n")


def abort_freeda():
    """Closes the GUI and aborts FREEDA pipeline"""

    root.quit()
    root.destroy()


def check_input():
    """Checks for essential user input"""

    ready = True

    if not wdirectory.get():
        logging.info("...FATAL_ERROR... : Choose working directory")
        ready = False

    if not clade.get():
        logging.info("...FATAL_ERROR... : Choose clade")
        ready = False

    if not threshold.get():
        logging.info("...FATAL_ERROR... : Choose search depth")
        ready = False

    all_proteins = [gene_name1.get(), gene_name2.get(), gene_name3.get()]
    if not any(all_proteins):
        logging.info("...FATAL_ERROR... : Choose at least one gene")
        ready = False

    return ready


def freeda_pipeline():
    """Main function running all freeda pipeline"""

    # use of the pipeline from GUI
    gui = True

    # deactivate the analyse button
    analyze_button["state"] = "disable"

    # get user input
    wdir = wdirectory.get() + "/"
    ref_species = clade.get()
    t = threshold.get()
    protein1 = gene_name1.get()
    protein2 = gene_name2.get()
    protein3 = gene_name3.get()

    # get user defined sites
    # gene 1
    gene1_site1_start = site11_start.get()
    gene1_site1_end = site11_end.get()
    gene1_site1_label = site11_label.get()

    gene1_site2_start = site12_start.get()
    gene1_site2_end = site12_end.get()
    gene1_site2_label = site12_label.get()

    gene1_site3_start = site13_start.get()
    gene1_site3_end = site13_end.get()
    gene1_site3_label = site13_label.get()

    # gene 2
    gene2_site1_start = site21_start.get()

    # gene 3
    gene3_site1_start = site31_start.get()

    global logging_window

    # change working directory to user indicated
    os.chdir(wdir)

    # ----------------------------------------#
    ######## GET USER INPUT ########
    # ----------------------------------------#

    # generate basic folders for input if not present
    input_extractor.generate_basic_folders(wdir)

    # get settings
    aligner = "mafft"

    # get all species and genome names
    all_genomes = [genome[1] for genome in genomes_preprocessing.get_names(ref_species, ref_genome=False)]

    # ----------------------------------------#
    ######## GET ALL INPUT DATA  ########
    # ----------------------------------------#

    # generate a reference Genome object
    ref_genome_present, ensembl, ref_species, ref_genomes_path, ref_genome_contigs_dict, \
        biotype, all_genes_ensembl = input_extractor.generate_ref_genome_object(wdir, ref_species)

    # check if provided gene names are present in ensembl object for ref assembly
    expected_proteins = [protein1, protein2, protein3]
    all_proteins = [protein for protein in expected_proteins if protein is not ""]

    if not input_extractor.validate_gene_names(all_proteins, all_genes_ensembl):
        return

    # stop pipeline if the reference genome is absent
    if not ref_genome_present:
        message = "\n...FATAL ERROR... : There is no reference genome detected -> exiting the pipeline now...\n"
        logging.info(message)
        #print("\n...FATAL ERROR... : There is no reference genome detected -> exiting the pipeline now...\n")
        return

    # get names of proteins

    for protein in all_proteins:

        message = "\n----------- * %s * -----------" % protein
        logging.info(message)
        #level = logging.WARNING
        #logger.log(level, message)   # DOESNT WORK
        #print("\n----------- * %s * -----------" % protein)
        # get structure prediction model from AlphaFold
        possible_uniprot_ids = input_extractor.get_uniprot_id(ref_species, protein)
        model_seq, uniprot_id = input_extractor.fetch_structure_prediction(wdir, ref_species,
                                                                               protein, possible_uniprot_ids)
        # get sequence input from ensembl
        input_correct, model_matches_input, microexon_present, microexons = input_extractor.extract_input(
            wdir, ref_species, ref_genomes_path,
            ref_genome_contigs_dict, ensembl, biotype,
            protein, model_seq, uniprot_id
        )

        if input_correct:
            message = "\nInput data have been generated for protein: %s\n\n" % protein
            logging.info(message)
            #print("\nInput data have been generated for protein: %s\n\n" % protein)

        if not input_correct:
            message = "\n...FATAL ERROR... : Input data generation FAILED for protein: %s " \
                      "- please remove from analysis -> exiting the pipeline now...\n" % protein
            logging.info(message)
            #print("\n...FATAL ERROR... : Input data generation FAILED for protein: %s - please remove from analysis"
            #      " -> exiting the pipeline now...\n" % protein)
            return

        if not model_matches_input:
            message = "\n...WARNING... : No matching structure prediction model is available for : %s " \
                  "-> cannot overlay FREEDA results onto a 3D structure\n" % protein
            logging.info(message)
            #print("...WARNING... : No matching structure prediction model is available for : %s "
            #      "-> cannot overlay FREEDA results onto a 3D structure\n" % protein)
            #print("...WARNING... : Protein will still be analyzed using PAML but without 3D structure overlay\n")

    # ----------------------------------------#
    ######## RUN BLAST ########
    # ----------------------------------------#

    print("Checking genome blast databases...")
    tblastn.run_blast(wdir, ref_species, all_proteins)

    # ----------------------------------------#
    ######## RUN EXON FINDING ########
    # ----------------------------------------#

    if exon_extractor.check_blast_output(wdir + "Blast_output/", t, all_proteins):
        result_path = exon_extractor.analyse_blast_results(wdir, wdir + "Blast_output/",
                                                           ref_species, int(t), all_proteins, all_genomes, aligner, gui,
                                                           logging_window)
    else:
        message = "\n     ...FATAL ERROR... : Genome of at least one species contains " \
              "no matches above the identity threshold used : %s -> use a lower one " \
              "-> exiting the pipeline now..." % t
        logging.info(message)
        #print("\n     ...FATAL ERROR... : Genome of at least one species contains "
        #      "no matches above the identity threshold used : %s -> use a lower one "
        #      "-> exiting the pipeline now..." % t)
        return

    # ----------------------------------------#
    ######## RUN PAML and PyMOL ########
    # ----------------------------------------#

    # run PAML
    nr_of_species_total_dict, PAML_logfile_name, day, failed_paml, \
            prots_under_pos_sel = paml_launcher.analyse_final_cds(wdir, ref_species, result_path,
                                                                  all_proteins, aligner, gui, logging_window)

    # visualize PAML result
    final_PAML_log_dict = paml_visualizer.analyse_PAML_results(wdir, result_path, all_proteins,
                                                               nr_of_species_total_dict, ref_species, PAML_logfile_name,
                                                               day, prots_under_pos_sel, failed_paml)
    get_results(final_PAML_log_dict)

    # run PyMOL
    for protein in all_proteins:

        # do not allow further analysis of failed paml runs
        if protein in failed_paml:
            continue

        # check if model seq and input seq match and check if exactly one model exists
        elif structure_builder.check_structure(wdir, ref_species, protein):
            successful = structure_builder.run_pymol(wdir, ref_species, result_path, protein, prots_under_pos_sel,
                                                                 offset=None)
            if not successful:
                message = "\nThe structure for : %s was not built successfully." % protein
                logging.info(message)
                #print("\nThe structure for : %s was not built successfully." % protein)
                continue
        else:
            message = "\nPrediction model for : %s DOES NOT match input sequence" \
                     "-> cannot overlay FREEDA results onto a 3D structure\n" % protein
            logging.info(message)
            #print("\nPrediction model for : %s DOES NOT match input sequence "
            #         "-> cannot overlay FREEDA results onto a 3D structure\n" % protein)

    logging.info("\nYou reached the end of FREEDA pipeline.")


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
    # accept only entry starting with one capital letter followed by small or big letters or numbers
    valid = re.match(r"^[1-9]{1}([0-9]{0,3}$)", residue) is not None  # max 4 digits, first one not a zero
    # button can be clicked only if gene names are valid
    analyze_button.state(["!disabled"] if valid else ["disabled"])
    # keystroke validation
    if op == "key":
        ok_so_far = re.match(r"^(?![\s\S])|[0-9]+$", residue) is not None
        if not ok_so_far:
            error_message2.set(message2)
        return ok_so_far
    elif op == "focusout":
        if not valid:
            error_message2.set(message2)
    return valid


def get_wdir():
    """Asks for a working directory."""
    directory = filedialog.askdirectory(initialdir=os.getcwd(), title="Select 'Data' folder")
    wdirectory.set(directory)


def get_results(final_PAML_log_dict):
    """Fetches key PAML results for each protein"""

    proteins = final_PAML_log_dict["Protein name"]
    lrts = final_PAML_log_dict["M8 vs M7 (LRT)"]
    pvalues = final_PAML_log_dict["M8 vs M7 (p-value)"]

    if gene_name1.get():
        g1_results_var.set(proteins.pop(0))
        g1_lrt_var.set(round(float(lrts.pop(0)), ndigits=3))
        pvalue = round(float(pvalues.pop(0)), ndigits=4)
        g1_pvalue_var.set(pvalue)
        if pvalue < 0.05:
            g1_pos_sel_var.set("YES")
            g1_pos_sel_entry.config(foreground="magenta")
        else:
            g1_pos_sel_var.set("NO")
            g1_pos_sel_entry.config(foreground="black")

    if gene_name2.get():
        g2_results_var.set(proteins.pop(0))
        g2_lrt_var.set(round(float(lrts.pop(0)), ndigits=3))
        pvalue = round(float(pvalues.pop(0)), ndigits=4)
        g2_pvalue_var.set(pvalue)
        if pvalue < 0.05:
            g2_pos_sel_var.set("YES")
            g2_pos_sel_entry.config(foreground="magenta")
        else:
            g2_pos_sel_var.set("NO")
            g2_pos_sel_entry.config(foreground="black")

    if gene_name3.get():
        g3_results_var.set(proteins.pop(0))
        g3_lrt_var.set(round(float(lrts.pop(0)), ndigits=3))
        pvalue = round(float(pvalues.pop(0)), ndigits=4)
        g3_pvalue_var.set(pvalue)
        if pvalue < 0.05:
            g3_pos_sel_var.set("YES")
            g3_pos_sel_entry.config(foreground="magenta")
        else:
            g3_pos_sel_var.set("NO")
            g3_pos_sel_entry.config(foreground="black")


def get_alignment():

    return


# set up the main window
root = Tk()
root.maxsize(width=1280, height=725)
#root.maxsize(width=100, height=100)
root.resizable(False, False)
#run = True

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
logo_frame.grid(column=0, row=0, columnspan=1, rowspan=3, sticky=(N, W, E, S), padx=5, pady=5)

# create settings frame
settings_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
settings_frame.grid(column=1, row=0, columnspan=3, rowspan=3, sticky=(N, W, E, S), padx=5, pady=5)

# create a gene 1 frame
gene1_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
gene1_frame.grid(column=0, row=11, columnspan=4, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
gene1_frame.columnconfigure(0, weight=3, minsize=50, uniform="group1")
gene1_frame.columnconfigure((1, 2), weight=1, minsize=50, uniform="group1")
gene1_frame.columnconfigure(3, weight=2, minsize=50, uniform="group1")

# create a gene 2 frame
gene2_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
gene2_frame.grid(column=0, row=14, columnspan=4, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
gene2_frame.columnconfigure(0, weight=3, uniform="group1")
gene2_frame.columnconfigure((1, 2), weight=1, uniform="group1")
gene2_frame.columnconfigure(3, weight=2, uniform="group1")

# create a gene 3 frame
gene3_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
gene3_frame.grid(column=0, row=17, columnspan=4, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
gene3_frame.columnconfigure(0, weight=3, uniform="group1")
gene3_frame.columnconfigure((1, 2), weight=1, uniform="group1")
gene3_frame.columnconfigure(3, weight=2, uniform="group1")

# create output frame
output_frame = ttk.Frame(mainframe, relief="sunken", padding="5 5 5 5")
output_frame.grid(column=1, row=0, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
output_frame.columnconfigure((0, 1, 2, 3), weight=1, uniform="group1")
output_frame.columnconfigure((4, 5, 6), weight=2, uniform="group1")

# create a logging window
logging_frame = ttk.Frame(output_frame, relief="ridge", padding="5 5 5 5")
logging_frame.grid(column=0, row=0, columnspan=7, sticky=(N, W, E, S), padx=5, pady=5)

# create results frame
results_labelframe = ttk.LabelFrame(output_frame, text="Results window")
results_labelframe.grid(column=0, row=1, columnspan=7, padx=5, pady=5, sticky=(N, W))
results_labelframe.columnconfigure((0, 1, 2, 3), weight=1, uniform="group1")
results_labelframe.columnconfigure((4, 5, 6), weight=2, uniform="group1")

# create results labels frame
results_labels = ttk.Frame(results_labelframe, padding="5 5 5 5")
results_labels.grid(column=0, row=3, columnspan=7, sticky=(N, W, E, S), padx=5, pady=5)
results_labels.columnconfigure((0, 1, 2, 3), weight=1, uniform="group1")
results_labels.columnconfigure((4, 5, 6), weight=2, uniform="group1")

# create user gene 1 result frame
g1_results_frame = ttk.Frame(results_labelframe, relief="ridge", padding="5 5 5 5")
g1_results_frame.grid(column=0, row=4, columnspan=8, sticky=(N, W, E, S), padx=5, pady=5)
g1_results_frame.columnconfigure((0, 1, 2, 3), weight=1, uniform="group1")
g1_results_frame.columnconfigure((4, 5, 6), weight=2, uniform="group1")

# create user gene 2 result frame
g2_results_frame = ttk.Frame(results_labelframe, relief="ridge", padding="5 5 5 5")
g2_results_frame.grid(column=0, row=5, columnspan=7, sticky=(N, W, E, S), padx=5, pady=5)
g2_results_frame.columnconfigure((0, 1, 2, 3), weight=1, uniform="group1")
g2_results_frame.columnconfigure((4, 5, 6), weight=2, uniform="group1")

# create user gene 3 result frame
g3_results_frame = ttk.Frame(results_labelframe, relief="ridge", padding="5 5 5 5")
g3_results_frame.grid(column=0, row=6, columnspan=7, sticky=(N, W, E, S), padx=5, pady=5)
g3_results_frame.columnconfigure((0, 1, 2, 3), weight=1, uniform="group1")
g3_results_frame.columnconfigure((4, 5, 6), weight=2, uniform="group1")


# CHECKS
check_gene_name_wrapper = (settings_frame.register(check_gene_name), "%P", "%V")
check_functional_residues = (settings_frame.register(check_functional_residues), "%P", "%V")

# ERRORS
error_message1 = StringVar()
message1 = "Invalid gene name (follow pattern: Cenpo \nfor rodents and : CENPO for others)"
error_message2 = StringVar()
message2 = "Invalid residue number (follow pattern: 100)"
# error labels
error_label1 = ttk.Label(input_frame, font="TkSmallCaptionFont", foreground="magenta", textvariable=error_message1)
error_label1.grid(column=0, row=22, columnspan=4, padx=5, pady=2, sticky="w")
error_label2 = ttk.Label(input_frame, font="TkSmallCaptionFont", foreground="magenta", textvariable=error_message2)
error_label2.grid(column=0, row=23, columnspan=4, padx=5, pady=2, sticky="w")

# LOGO
canvas = Canvas(logo_frame, width=200, height=200)
canvas.pack()
freeda_img = ImageTk.PhotoImage(Image.open(os.getcwd() + "/freeda_img.png"))
canvas.create_image(100, 100, image=freeda_img)

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

# CLADE
clade = StringVar()
ttk.Label(settings_frame, text="Clade").grid(column=0, row=0, pady=5, sticky=(W))
rodents = ttk.Radiobutton(settings_frame, text="Rodents (mouse <-> rat)", variable=clade, value="Mm")
rodents.grid(column=0, row=1, columnspan=2, sticky=(W))
primates = ttk.Radiobutton(settings_frame, text="Primates (human <-> spider monkey)", variable=clade, value="Hs")
primates.grid(column=0, row=2, columnspan=2, sticky=(W))
carnivores = ttk.Radiobutton(settings_frame, text="Carnivores (cat <-> dog)", variable=clade, value="Cf")
carnivores.grid(column=0, row=3, columnspan=2, sticky=(W))
birds = ttk.Radiobutton(settings_frame, text="Birds (chicken <-> quail)", variable=clade, value="Gg")
birds.grid(column=0, row=4, columnspan=2, sticky=(W))

# BLAST THRESHOLD
threshold = IntVar()
ttk.Label(settings_frame, text="Blast search").grid(column=0, row=5, pady=5, sticky=(W))
ttk.Radiobutton(settings_frame, text="Deep", variable=threshold,
                       value=30).grid(column=0, row=6, sticky=(W))
ttk.Radiobutton(settings_frame, text="Medium (recommended)", variable=threshold,
                         value=50).grid(column=0, row=7, sticky=(W))
ttk.Radiobutton(settings_frame, text="Shallow", variable=threshold,
                          value=70).grid(column=0, row=8, sticky=(W))

# RESIDUES
ttk.Label(input_frame, text="Indicate functional residues -----> ").grid(column=0, row=9, columnspan=2,
                                                                     sticky=(W, E), padx=5)
ttk.Label(input_frame, text="start").grid(column=1, row=9)
ttk.Label(input_frame, text="end").grid(column=2, row=9)
ttk.Label(input_frame, text="label").grid(column=3, row=9)

# USER INPUT GENE 1
gene_name1 = StringVar()
ttk.Label(gene1_frame, text="Gene name").grid(column=0, row=11, padx=6, sticky=(W))
name1 = ttk.Entry(gene1_frame, textvariable=gene_name1, validate="all", validatecommand=check_gene_name_wrapper)
name1.grid(column=0, row=12, padx=5, pady=2, sticky=(W))
dup1 = StringVar()
ttk.Checkbutton(gene1_frame, text="Duplication expected", #command=duplication_expected,
                              variable=dup1, onvalue="on", offvalue="off").grid(column=0, row=13, padx=6, sticky=(W))

s11_start = StringVar()
s11_end = StringVar()
s11_label = StringVar()
site11_start = ttk.Entry(gene1_frame, textvariable=s11_start, validate="all", validatecommand=check_functional_residues)
site11_start.grid(column=1, row=11, padx=5, pady=2, sticky=(W))
site11_end = ttk.Entry(gene1_frame, textvariable=s11_end, validate="all", validatecommand=check_functional_residues)
site11_end.grid(column=2, row=11, padx=5, pady=2, sticky=(W))
site11_label = ttk.Entry(gene1_frame, textvariable=s11_label, validate="all")
site11_label.grid(column=3, row=11, padx=5, pady=2, sticky=(W))

s12_start = StringVar()
s12_end = StringVar()
s12_label = StringVar()
site12_start = ttk.Entry(gene1_frame, textvariable=s12_start, validate="all", validatecommand=check_functional_residues)
site12_start.grid(column=1, row=12, padx=5, pady=2, sticky=(W))
site12_end = ttk.Entry(gene1_frame, textvariable=s12_end, validate="all", validatecommand=check_functional_residues)
site12_end.grid(column=2, row=12, padx=5, pady=2, sticky=(W))
site12_label = ttk.Entry(gene1_frame, textvariable=s12_label, validate="all")
site12_label.grid(column=3, row=12, padx=5, pady=2, sticky=(W))

s13_start = StringVar()
s13_end = StringVar()
s13_label = StringVar()
site13_start = ttk.Entry(gene1_frame, textvariable=s13_start, validate="all", validatecommand=check_functional_residues)
site13_start.grid(column=1, row=13, padx=5, pady=2, sticky=(W))
site13_end = ttk.Entry(gene1_frame, textvariable=s13_end, validate="all", validatecommand=check_functional_residues)
site13_end.grid(column=2, row=13, padx=5, pady=2, sticky=(W))
site13_label = ttk.Entry(gene1_frame, textvariable=s13_label, validate="all")
site13_label.grid(column=3, row=13, padx=5, pady=2, sticky=(W))

# USER INPUT GENE 2
gene_name2 = StringVar()
ttk.Label(gene2_frame, text="Gene name").grid(column=0, row=14, padx=6, sticky=(W))
name2 = ttk.Entry(gene2_frame, textvariable=gene_name2, validate="all", validatecommand=check_gene_name_wrapper)
name2.grid(column=0, row=15, padx=5, pady=5, sticky=(W))
dup2 = StringVar()
ttk.Checkbutton(gene2_frame, text="Duplication expected", #command=duplication_expected,
                              variable=dup2, onvalue="on", offvalue="off").grid(column=0, row=16, padx=6, sticky=(W))

s21_start = StringVar()
s21_end = StringVar()
s21_label = StringVar()
site21_start = ttk.Entry(gene2_frame, textvariable=s21_start, validate="all", validatecommand=check_functional_residues)
site21_start.grid(column=1, row=14, padx=5, pady=2, sticky=(W))
site21_end = ttk.Entry(gene2_frame, textvariable=s21_end, validate="all", validatecommand=check_functional_residues)
site21_end.grid(column=2, row=14, padx=5, pady=2, sticky=(W))
site21_label = ttk.Entry(gene2_frame, textvariable=s21_label, validate="all")
site21_label.grid(column=3, row=14, padx=5, pady=2, sticky=(W))

s22_start = StringVar()
s22_end = StringVar()
s22_label = StringVar()
site22_start = ttk.Entry(gene2_frame, textvariable=s22_start, validate="all", validatecommand=check_functional_residues)
site22_start.grid(column=1, row=15, padx=5, pady=2, sticky=(W))
site22_end = ttk.Entry(gene2_frame, textvariable=s22_end, validate="all", validatecommand=check_functional_residues)
site22_end.grid(column=2, row=15, padx=5, pady=2, sticky=(W))
site22_label = ttk.Entry(gene2_frame, textvariable=s22_label, validate="all")
site22_label.grid(column=3, row=15, padx=5, pady=2, sticky=(W))

s23_start = StringVar()
s23_end = StringVar()
s23_label = StringVar()
site23_start = ttk.Entry(gene2_frame, textvariable=s23_start, validate="all", validatecommand=check_functional_residues)
site23_start.grid(column=1, row=16, padx=5, pady=2, sticky=(W))
site23_end = ttk.Entry(gene2_frame, textvariable=s23_end, validate="all", validatecommand=check_functional_residues)
site23_end.grid(column=2, row=16, padx=5, pady=2, sticky=(W))
site23_label = ttk.Entry(gene2_frame, textvariable=s23_label, validate="all")
site23_label.grid(column=3, row=16, padx=5, pady=2, sticky=(W))

# USER INPUT GENE 3
gene_name3 = StringVar()
ttk.Label(gene3_frame, text="Gene name").grid(column=0, row=17, padx=6, sticky=(W))
name3 = ttk.Entry(gene3_frame, textvariable=gene_name3, validate="all", validatecommand=check_gene_name_wrapper)
name3.grid(column=0, row=18, padx=5, pady=5, sticky=(W))
dup3 = StringVar()
ttk.Checkbutton(gene3_frame, text="Duplication expected", #command=duplication_expected,
                              variable=dup3, onvalue="on", offvalue="off").grid(column=0, row=19, padx=6, sticky=(W))

s31_start = StringVar()
s31_end = StringVar()
s31_label = StringVar()
site31_start = ttk.Entry(gene3_frame, textvariable=s31_start, validate="all", validatecommand=check_functional_residues)
site31_start.grid(column=1, row=17, padx=5, pady=2, sticky=(W))
site31_end = ttk.Entry(gene3_frame, textvariable=s31_end, validate="all", validatecommand=check_functional_residues)
site31_end.grid(column=2, row=17, padx=5, pady=2, sticky=(W))
site31_label = ttk.Entry(gene3_frame, textvariable=s31_label, validate="all")
site31_label.grid(column=3, row=17, padx=5, pady=2, sticky=(W))

s32_start = StringVar()
s32_end = StringVar()
s32_label = StringVar()
site32_start = ttk.Entry(gene3_frame, textvariable=s32_start, validate="all", validatecommand=check_functional_residues)
site32_start.grid(column=1, row=18, padx=5, pady=2, sticky=(W))
site32_end = ttk.Entry(gene3_frame, textvariable=s32_end, validate="all", validatecommand=check_functional_residues)
site32_end.grid(column=2, row=18, padx=5, pady=2, sticky=(W))
site32_label = ttk.Entry(gene3_frame, textvariable=s32_label, validate="all")
site32_label.grid(column=3, row=18, padx=5, pady=2, sticky=(W))

s33_start = StringVar()
s33_end = StringVar()
s33_label = StringVar()
site33_start = ttk.Entry(gene3_frame, textvariable=s33_start, validate="all", validatecommand=check_functional_residues)
site33_start.grid(column=1, row=19, padx=5, pady=2, sticky=(W))
site33_end = ttk.Entry(gene3_frame, textvariable=s33_end, validate="all", validatecommand=check_functional_residues)
site33_end.grid(column=2, row=19, padx=5, pady=2, sticky=(W))
site33_label = ttk.Entry(gene3_frame, textvariable=s33_label, validate="all")
site33_label.grid(column=3, row=19, padx=5, pady=2, sticky=(W))

# BUTTONS

# Working directory button
wdir_button = ttk.Button(input_frame, text="Set working directory", command=get_wdir)
wdir_button.grid(column=0, row=21, sticky=(N, W, E, S), padx=5, pady=2)
wdirectory = StringVar()
wdir_entry = ttk.Entry(input_frame, textvariable=wdirectory)
wdir_entry.grid(column=1, row=21, columnspan=4, sticky=(N, W, E, S), padx=5, pady=2)

# ANALYZE button
analyze_button = ttk.Button(input_frame, text="Analyze", state="normal", command=thread_freeda)  # default="active"
analyze_button.grid(column=1, row=22, columnspan=2, padx=5, pady=2, sticky=(E))

# ABORT button
abort_button = ttk.Button(input_frame, text="ABORT", state="normal", command=abort_freeda)  # default="active"
abort_button.grid(column=3, row=22, padx=5, pady=2, sticky=(W))

# RESULTS button
results_button = ttk.Button(results_labelframe, text="Results sheet", state="normal",
                            command=lambda x="Results sheet": get_results(x))
results_button.grid(column=0, row=7, columnspan=3, sticky=(W), padx=5, pady=5)

# ALIGNMENT BUTTONS
alignment_button1 = ttk.Button(g1_results_frame, text="Alignment", state="normal", command=get_results)
alignment_button1.grid(column=4, row=0, sticky=(W, E), padx=5, pady=5)
alignment_button2 = ttk.Button(g2_results_frame, text="Alignment", state="normal", command=get_results)
alignment_button2.grid(column=4, row=0, sticky=(W, E), padx=5, pady=5)
alignment_button3 = ttk.Button(g3_results_frame, text="Alignment", state="normal", command=get_results)
alignment_button3.grid(column=4, row=0, sticky=(W, E), padx=5, pady=5)

# PAML graph buttons
PAML_graph_button1 = ttk.Button(g1_results_frame, text="PAML graph", state="normal", command=get_results)
PAML_graph_button1.grid(column=5, row=0, sticky=(W, E), padx=5, pady=5)
PAML_graph_button2 = ttk.Button(g2_results_frame, text="PAML graph", state="normal", command=get_results)
PAML_graph_button2.grid(column=5, row=0, sticky=(W, E), padx=5, pady=5)
PAML_graph_button3 = ttk.Button(g3_results_frame, text="PAML graph", state="normal", command=get_results)
PAML_graph_button3.grid(column=5, row=0, sticky=(W, E), padx=5, pady=5)

# Structure buttons
structure_button1 = ttk.Button(g1_results_frame, text="Structure", state="normal", command=get_results)
structure_button1.grid(column=6, row=0, sticky=(W, E), padx=5, pady=5)
structure_button2 = ttk.Button(g2_results_frame, text="Structure", state="normal", command=get_results)
structure_button2.grid(column=6, row=0, sticky=(W, E), padx=5, pady=5)
structure_button3 = ttk.Button(g3_results_frame, text="Structure", state="normal", command=get_results)
structure_button3.grid(column=6, row=0, sticky=(W, E), padx=5, pady=5)


# LABELS

#ttk.Label(results_frame, text="Results window").grid(column=0, row=0, sticky=(W))
ttk.Label(results_labels, text="Gene").grid(column=0, row=3, padx=5)
ttk.Label(results_labels, text="Pos. sel.").grid(column=1, row=3, padx=5)
ttk.Label(results_labels, text="LRT").grid(column=2, row=3, padx=5)
ttk.Label(results_labels, text="p-value").grid(column=3, row=3, padx=5)

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




#ttk.Label(mainframe, text="Gene name").grid(column=0, row=11, sticky=(W))

#threshold = StringVar()
#threshold_label = ttk.Label(mainframe, textvariable=threshold).grid(column=1, row=1, sticky=(E))
#threshold_entry = ttk.Entry(mainframe, textvariable=threshold, width=2).grid(column=1, row=1, sticky=(W))


# make button to run the pipeline
#button = ttk.Button(mainframe, text="Analyze", default="active", command=freeda.freeda_pipeline)
# will execute script attached to the button when left mouse clicked
#root.bind("<ButtonPress-1>", lambda e: button.invoke())


#root.after(1000, freeda_pipeline)

root.mainloop()
