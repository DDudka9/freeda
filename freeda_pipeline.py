#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copyright 2022 - Damian Dudka and R. Brian Akins - contact: damiandudka0@gmail.com

This file is part of FREEDA.

FREEDA is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

FREEDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with FREEDA.
If not, see <https://www.gnu.org/licenses/>.

"""


"""

Main module of the FREEDA package used for debugging. Use from terminal only.

"""

print("\nImporting all modules and libraries...\n")

from freeda import folder_generator
from freeda import input_extractor
from freeda import tblastn
from freeda import exon_extractor
from freeda import paml_launcher
from freeda import paml_visualizer
from freeda import structure_builder
from freeda import genomes_preprocessing
import os
import shutil
from sys import platform
import logging


def freeda_pipeline(wdir=None, ref_species=None, t=None, codon_frequencies=None, excluded_species=None):
    """Main function running all freeda pipeline from command line"""

    if wdir is None:
        wdir = os.getcwd() + "/"

    if ref_species is None:
        ref_species = "Mm"
    elif ref_species not in ["Mm", "Hs", "Rn", "Cf", "Gg", "Fc", "Dme"]:
        print("\n...FATAL_ERROR... : Invalid reference species (try: Mm for mouse, Hs for human) "
              "-> exiting the pipeline now...")
        return

    # check if pymol is installed
    if not shutil.which("pymol"):
        print("\n...FATAL_ERROR... : Pymol not installed "
                     "\n macOS -> follow README file or go to https://pymol.org/2/ to download and install Pymol"
                     "\n ubuntu -> follow README file or go to Software Manager and download and install Pymol")
        return

    # initial percent identity threshold for blast matches analysis
    if t is None:
        t = 60

    if codon_frequencies is None:
        codon_frequencies = "F3X4"
    elif codon_frequencies.upper() != "F3X4" and codon_frequencies.upper() != "F61" \
                                            and codon_frequencies.upper() != "F3X4, F61":
        print("\n...FATAL_ERROR... : Invalid codon frequencies (try: F3X4 or F61 or F3X4, F61) "
              "-> exiting the pipeline now...")
        return

    if excluded_species is None:
        final_excluded_species = {}

    if excluded_species is not None:
        excluded_species = str(excluded_species).split(" ")
        final_excluded_species = {}
        available_species = genomes_preprocessing.get_available_species(ref_species)
        for ex_species in excluded_species:
            # check if users input matches available species
            if ex_species in available_species:
                final_excluded_species[ex_species] = available_species[ex_species]  # add genome as value to species key

    user_input0 = None
    user_input1 = None
    user_input2 = None
    user_input3 = None

    # ----------------------------------------#
    ######## GET USER INPUT ########
    # ----------------------------------------#

    print("Choose which parts of the pipeline you would like to run (all 'y' is a good strategy for single genes) : ")
    while user_input0 != "y" and user_input0 != "n":
        user_input0 = input("\n(FREEDA) Should I get input data automatically? (y / n)\n").lower()
        if user_input0.lower() != "y" and user_input0.lower() != "n":
            print("Please answer y or n\n")

    while user_input1 != "y" and user_input1 != "n":
        user_input1 = input("(FREEDA) Should I run blast with the input data? (y / n)\n").lower()
        if user_input1.lower() != "y" and user_input1.lower() != "n":
            print("Please answer y or n\n")

    while user_input2 != "y" and user_input2 != "n":
        user_input2 = input("(FREEDA) Should I find exons based on the blast results? (y / n)\n").lower()
        if user_input2.lower() != "y" and user_input2.lower() != "n":
            print("Please answer y or n\n")

    while user_input3 != "y" and user_input3 != "n":
        user_input3 = input("(FREEDA) Should I perform molecular evolution analysis (PAML)? (y / n)\n").lower()
        if user_input3.lower() != "y" and user_input3.lower() != "n":
            print("Please answer y or n\n")

    # generate basic folders for input if not present
    folder_generator.generate_basic_folders(wdir)
    # get all genes to be analyzed (these are gene names)
    all_genes = [gene.rstrip("\n") for gene in open(wdir + "genes.txt", "r").readlines() if gene != "\n"]

    # get all species and genome names
    all_genomes = [genome[1] for genome in genomes_preprocessing.get_names(wdir, ref_species, final_excluded_species)]

    # check if the user had previously obtained data for given list of genes
    if user_input0 == "n":
        input_present = True

        for gene in all_genes:

            structure_path = wdir + "Structures/" + gene + "_" + ref_species
            if "model_matches_input_seq.txt" in os.listdir(structure_path) \
                    or "model_incompatible.txt" in os.listdir(structure_path):
                print("\nAll input data and structure model for : %s are present." % gene)
            else:
                input_present = False
                print("\n...WARNING... : Data are missing for : %s" % gene)
        if not input_present:
            print("\n...FATAL_ERROR... : Input data for at least one gene are missing "
                  "-> exiting the pipeline now...")
            return

    # ----------------------------------------#
    ######## CONDITIONS NOT ALLOWED ########
    # ----------------------------------------#

    # forgot about blast
    if user_input0 == "y" and user_input1 == "n" and user_input2 == "y" and user_input3 == "n":
        print("\n...FATAL ERROR... : You need to perform blast before exon finding.")
        return

    # forgot about blast or exon finding
    if user_input0 == "y" and (user_input1 == "n" or user_input2 == "n") and user_input3 == "y":
        print("\n...FATAL ERROR... : You need to perform blast and exon finding before molecular evolution analysis.")
        return

    # forgot about exon finding
    if user_input0 == "n" and user_input1 == "y" and user_input2 == "n" and user_input3 == "y":
        print("\n...FATAL ERROR... : You need to perform exon finding before molecular evolution analysis.")
        return

    # ----------------------------------------#
    ######## INSTALL PYMOL IF NEEDED ########
    # ----------------------------------------#

    if not shutil.which("pymol"):
        if platform == "linux" or platform == "linux2":
            logging.info("\nPyMOL not found in the PATH.")
            pymol_choice = input("\nInstall PyMOL into current working directory? (y/N)")
            if pymol_choice == "y":
                structure_builder.install_pymol_linux(wdir)
            else:
                logging.info("\nPyMOL is required to run FREEDA. It can be installed using a package manager. "
                             "For example, try:\n    sudo apt-get install pymol")
                return
        elif platform == "darwin":
            message = "\nPyMOL not found in the path. Install PyMOL using homebrew with:\n    brew install pymol"
            logging.info(message)
            return

    # ----------------------------------------#
    ######## GET ALL INPUT DATA  ########
    # ----------------------------------------#

    # User wants to generate input data
    if user_input0 == "y":

        # generate a reference Genome object
        # WHY DOES IT RETURN REF_SPECIES?
        ref_genome_present, ensembl, ref_species, ref_genomes_path, ref_genome_contigs_dict, \
                        biotype, all_genes_ensembl, all_gene_ids_ensembl = input_extractor.generate_ref_genome_object(wdir, ref_species)

        # check if provided gene names are present in ensembl object for ref assembly
        if not input_extractor.validate_gene_names(ref_species, all_genes, all_genes_ensembl, ensembl):
            return

        # stop pipeline if the reference genome is absent
        if not ref_genome_present:
            print("\n...FATAL ERROR... : There is no reference genome detected -> exiting the pipeline now...")
            return

        # get names of

        for gene in all_genes:

            print("\n----------- * %s * -----------" % gene)
            # get structure prediction model from AlphaFold
            possible_uniprot_ids = input_extractor.get_uniprot_id(ref_species, gene, ensembl)
            model_seq, uniprot_id = input_extractor.fetch_structure_prediction(wdir, ref_species,
                                                                               gene, possible_uniprot_ids)
            # get sequence input from ensembl
            input_correct, model_matches_input, microexon_present, microexons = input_extractor.extract_input(
                wdir, ref_species, ref_genomes_path,
                ref_genome_contigs_dict, ensembl, biotype,
                gene, model_seq, uniprot_id
            )

            if input_correct:
                print("\nInput data have been generated for gene: %s\n\n" % gene)

            if not input_correct:
                print("\n...FATAL ERROR... : Input data generation FAILED for gene: %s - please remove %s from analysis"
                      " -> exiting the pipeline now...\n" % (gene, gene))
                return

            if not model_matches_input:
                print("...WARNING... : No entry in searched Ensembl release matches prediction model for : %s [%s] "
                      "-> cannot overlay FREEDA results onto a 3D structure\n" % (gene, uniprot_id))
                print("...WARNING... : Gene will still be analyzed using PAML but without 3D structure overlay\n")


    # ----------------------------------------#
    ######## RUN BLAST ########
    # ----------------------------------------#

    blast_output_path = wdir + "Blast_output/"
    if user_input1 == "y":
        print("Checking genome blast databases...")
        tblastn.run_blast(wdir, ref_species, all_genes, final_excluded_species)

    # ----------------------------------------#
    ######## RUN EXON FINDING ########
    # ----------------------------------------#

    if user_input2 == "y":
        if exon_extractor.check_blast_output(ref_species, blast_output_path, t, all_genes):
            result_path = exon_extractor.analyze_blast_results(wdir, blast_output_path, ref_species, int(t),
                                                               all_genes, all_genomes,
                                                               final_excluded_species)
        else:
            return

    # ----------------------------------------#
    ######## RUN PAML and PyMOL ########
    # ----------------------------------------#

    if user_input3 == "y" and user_input2 == "n":

        nr_of_tries = 1
        while nr_of_tries <= 3:
            user_input4 = input("----->  Indicate folder with results (no slashes or quotes): ")
            result_path = wdir + user_input4 + "/Raw_data/"

            if os.path.isdir(result_path) is False:
                print("\n(FREEDA) I could not find your results folder (%s/3)" % nr_of_tries)
                nr_of_tries += 1
            else:
                nr_of_tries = float("inf")
                result_path = wdir + user_input4 + "/Raw_data/"

                # run PAML
                nr_of_species_total_dict, PAML_logfile_name, day, \
                failed_paml, genes_under_pos_sel = paml_launcher.analyze_final_cds(wdir, ref_species,
                                                            result_path, all_genes, codon_frequencies)

                # visualize PAML result
                paml_visualizer.analyze_PAML_results(wdir, result_path, all_genes,
                                                     nr_of_species_total_dict, ref_species, PAML_logfile_name,
                                                     day, genes_under_pos_sel, failed_paml, codon_frequencies)
                # run PyMOL
                for gene in all_genes:

                    # do not allow further analysis of failed paml runs
                    if gene in failed_paml:
                        continue

                    # check if model seq and input seq match and check if exactly one model exists
                    elif structure_builder.check_structure(wdir, ref_species, gene):
                        successful = structure_builder.run_pymol(wdir, ref_species, result_path,
                                                                 gene, genes_under_pos_sel)
                        if not successful:
                            print("\nThe structure for : %s was not built successfully." % gene)
                            continue
                    else:
                        print("\nPrediction model for : %s DOES NOT match input sequence "
                              "-> cannot overlay FREEDA results onto a 3D structure\n" % gene)

    if user_input3 == "y" and user_input2 == "y":

        # run PAML
        nr_of_species_total_dict, PAML_logfile_name, day, \
        failed_paml, genes_under_pos_sel = paml_launcher.analyze_final_cds(wdir, ref_species,
                                                            result_path, all_genes, codon_frequencies)

        # visualize PAML result
        paml_visualizer.analyze_PAML_results(wdir, result_path, all_genes,
                                             nr_of_species_total_dict, ref_species, PAML_logfile_name,
                                             day, genes_under_pos_sel, failed_paml, codon_frequencies)
        # run PyMOL
        for gene in all_genes:

            # do not allow further analysis of failed paml runs
            if gene in failed_paml:
                continue

            # check if model seq and input seq match and check if exactly one model exists
            elif structure_builder.check_structure(wdir, ref_species, gene):
                successful = structure_builder.run_pymol(wdir, ref_species, result_path,
                                                         gene, genes_under_pos_sel)
                if not successful:
                    print("\nThe structure for : %s was not built successfully." % gene)
                    continue

            else:
                print("\nPrediction model for : %s DOES NOT match input sequence "
                      "-> cannot overlay FREEDA results onto a 3D structure\n" % gene)

    print("\nYou reached the end of FREEDA pipeline.")


# ----------------------------------------#
######## RUN as command line ########
# ----------------------------------------#

# need to provide an absolute path to the main when running in command line:
# (py37) python /Users/damian/PycharmProjects/freeda_2.0/freeda_pipeline.py -d /Volumes/DamianEx_2/Data/ -rs "Mm" -t 30
# you can get abs path using:
# (base) brew install coreutils
# (base) realpath freeda_pipeline.py

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-d", "--wdir",
                        help="specify working directory (absolute path to Data folder ex. /Users/user/Data/)", type=str,
                        default=None)
    parser.add_argument("-rs", "--ref_species",
                        help="specify reference organism (default is mouse)", type=str, default="Hs")
    parser.add_argument("-t", "--blast_threshold",
                        help="specify percentage identity threshold for blast (default is 60)", type=int, default=30)
    parser.add_argument("-f", "--codon_frequencies",
                        help="specify codon frequency models (F3x4 is default)", type=str, default="F3X4")
    parser.add_argument("-es", "--excluded_species",
                        help="specify species to exclude (e.g. Ha Gs)", type=str, default="")

    args = parser.parse_args()
    freeda_pipeline(wdir=args.wdir, ref_species=args.ref_species, t=args.blast_threshold,
                    codon_frequencies=args.codon_frequencies, excluded_species=args.excluded_species)
