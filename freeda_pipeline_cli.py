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
from freeda import pyinstaller_compatibility
import os
import sys
import re
import shutil


def freeda_pipeline(wdir=None, ref_species=None, t=None, codon_frequencies=None, subgroup=None, excluded_species=None):
    """Main function running all freeda pipeline from command line"""

    gui = None
    logging_window = None

    # check if pymol is installed
    if not shutil.which("pymol"):
        print("\n...FATAL_ERROR... : Pymol not installed "
                     "\n macOS -> follow README file or go to https://pymol.org/2/ to download and install Pymol"
                     "\n ubuntu -> follow README file or go to Software Manager and download and install Pymol")
        return

    # ----------------------------------------#
    ######## ASSIGN ENVIRONMENT VARIABLES ########
    # ----------------------------------------#

    if pyinstaller_compatibility.is_bundled():
        os.environ["MAFFT_BINARIES"] = pyinstaller_compatibility.resource_path("mafft_bin")
        os.environ["REQUESTS_CA_BUNDLE"] = pyinstaller_compatibility.resource_path("certifi/cacert.pem")

    # ----------------------------------------#
    ######## INSTALL PYMOL IF NEEDED ########
    # ----------------------------------------#

    if not shutil.which("pymol") and \
            not os.path.exists(os.path.join("/", "Applications", "PyMOL.app", "Contents", "MacOS", "PyMOL")):
        if sys.platform == "linux" or sys.platform == "linux2":
            print("\nPyMOL not found in the PATH. Checking for PyMOL in the current working directory.")
            structure_builder.install_pymol_linux(wdir)
        else:
            message = "\nWELCOME! It seems that you are running FREEDA\nfor the first time -> " \
                      "Please download PyMOL from:\n\nhttps://pymol.org/\n\nand place it in the Applications folder\n" \
                      "Then close the FREEDA application (or click ABORT) and open it again"
            print(message)
            return

    # ----------------------------------------#
    ######## GET USER INPUT ########
    # ----------------------------------------#

    # SET WORKING DIRECTORY
    if wdir is None:
        wdir = input("\n(FREEDA) Pick working directory (/Users/user/Data/):\n")
        if not os.path.exists(wdir):
            print("\nYou did not pick a valid directory -> exiting the pipeline...\n")
            return
        else:
            os.chdir(wdir)

    # SET REFERENCE SPECIES
    if ref_species is None:
        ref_species_asw = input("\n(FREEDA) Pick reference species (mouse, human, dog, chicken, fly):\n").lower()
        if ref_species_asw not in ["mouse", "human", "dog", "chicken", "fly"]:
            print("\nYou did not pick a valid reference species -> exiting the pipeline...\n")
            return
        else:
            if ref_species_asw == "mouse":
                ref_species = "Mm"
            elif ref_species_asw == "human":
                ref_species = "Hs"
            elif ref_species_asw == "dog":
                ref_species = "Cf"
            elif ref_species_asw == "chicken":
                ref_species = "Gg"
            else:
                ref_species = "Dme"

    # SET SUBGROUP
    if subgroup is None:
        subgroup = input("\n(FREEDA) (optional) Pick subgroup (hominoidea, catarrhini, caniformia, melanogaster):\n").lower()
        if subgroup not in ["hominoidea", "catarrhini", "caniformia", "melanogaster"]:
            print("You did not pick a valid subgroup (optional)\n")
            subgroup = ""

    # SET EXCLUDED SPECIES
    if excluded_species is None:
        excluded_species = input("\n(FREEDA) (optional) Pick species to exclude (e.g. Ha Gd):\n")
        excluded_species = str(excluded_species).split(" ")
        final_excluded_species = {}
        available_species = genomes_preprocessing.get_available_species(ref_species)
        for ex_species in excluded_species:
            # check if users input matches available species
            if ex_species in available_species:
                final_excluded_species[ex_species] = available_species[ex_species]  # add genome as value to species key
    else:
        print("You did not pick any species to exclude\n")
        final_excluded_species = {}


    # SET CODON FREQUENCIES
    if codon_frequencies is None:
        codon_frequencies_ans = input("\n(FREEDA) (optional) Pick a number to set codon frequencies "
                                      "(1 - F3X4, "
                                      "2 - F61, "
                                      "3 - F3X4 and F61):\n")
        if codon_frequencies_ans == "1":
            codon_frequencies = "F3X4"
        elif codon_frequencies_ans == "2":
            codon_frequencies = "F61"
        elif codon_frequencies_ans == "3":
            codon_frequencies = "F3X4, F61"
        else:
            print("You did not pick a valid number (F3X4 remains a default)\n")
            codon_frequencies = "F3X4"

    # SET GENE PARAMETERS
    all_genes_dict = {}
    user_gene = ""
    advanced_option_1 = "Advanced options (OFF)"
    advanced_option_2 = "Advanced options (OFF)"
    site1_label = ""
    site2_label = ""
    site3_label = ""
    site1_start = ""
    site2_start = ""
    site3_start = ""
    site1_end = ""
    site2_end = ""
    site3_end = ""

    # gene parameters default
    all_genes_dict[user_gene] = [advanced_option_1, advanced_option_2,
                                 [site1_label, site1_start, site1_end],
                                 [site2_label, site2_start, site2_end],
                                 [site3_label, site3_start, site3_end]]

    user_genes = input("\n(FREEDA) Pick genes to analyze (e.g. CENPO MX1):\n").split(" ")

    for user_gene in user_genes:

        if not re.match(r"^[A-Z]{1}([A-Za-z0-9]+$)", user_gene):
            print("\nYou did not pick any genes or used invalid characters -> exiting the pipeline...\n")
            return

        else:
            # define advanced options and labels
            advanced_options = input("\n(FREEDA) (optional) Access advanced options for gene: %s "
                                         "(yes or no)\n" % user_gene).lower()
            if advanced_options == "yes":
                advanced_option_1 = input("\n(FREEDA) (advanced) For gene: %s pick one number: "
                                               "1 - Duplication expected, "
                                               "2 - Tandem duplication expected, "
                                               "3 - Long introns expected (>50kb)\n" % user_gene)
                if advanced_option_1 == "1":
                    advanced_option_1 = "Duplication expected"
                elif advanced_option_1 == "2":
                    advanced_option_1 = "Tandem duplication expected"
                elif advanced_option_1 == "3":
                    advanced_option_1 = "Long introns expected (>50kb)"
                else:
                    print("\nYou did not select any of the available options for gene: %s\n" % user_gene)
                    advanced_option_1 = "Advanced options (OFF)"

                advanced_option_2 = input("\n(FREEDA) (advanced) Do you expect common domains in gene: %s? "
                                          "(yes or no)\n"
                                              % user_gene).lower()
                if advanced_option_2 == "yes":
                    advanced_option_2 = "Common domains expected"
                else:
                    print("\nYou do not expect common domains in gene: %s\n" % user_gene)
                    advanced_option_2 = "Advanced options (OFF)"

            # no advanced options selected needed for this gene
            else:
                print("\nYou did not not select any advanced options for gene: %s\n" % user_gene)
                # gene advanced options default
                all_genes_dict[user_gene] = [advanced_option_1, advanced_option_2,
                                             [site1_label, site1_start, site1_end],
                                             [site2_label, site2_start, site2_end],
                                             [site3_label, site3_start, site3_end]]

            # define labels
            site_label = input("\n(FREEDA) (optional) Would you like to label up to 3 regions in gene: %s "
                                       "(yes or no)\n" % user_gene).lower()
            if site_label == "yes":
                site1_label = input("\n(FREEDA) (optional) Label first region\n")
                if site1_label != "":
                    site1_start = input("\n(FREEDA) (optional) Pick first residue number for region: %s\n" % site1_label)
                    if not re.match(r"^([0-9]+$)", site1_start):
                        print("\nInvalid residue number for region: %s\n" % site1_label)
                    site1_end = input("\n(FREEDA) (optional) Pick last residue number for region: %s\n" % site1_label)
                    if not re.match(r"^([0-9]+$)", site1_end):
                        print("\nInvalid residue number for region: %s\n" % site1_label)

                site2_label = input("\n(FREEDA) (optional) Label second region\n")
                if site2_label != "":
                    site2_start = input("\n(FREEDA) (optional) Pick first residue number for region: %s\n" % site2_label)
                    if not re.match(r"^([0-9]+$)", site2_start):
                        print("\nInvalid residue number for region: %s\n" % site2_label)
                    site2_end = input("\n(FREEDA) (optional) Pick last residue number for region: %s\n" % site2_label)
                    if not re.match(r"^([0-9]+$)", site2_end):
                        print("\nInvalid residue number for region: %s\n" % site2_label)

                site3_label = input("\n(FREEDA) (optional) Label third region\n")
                if site3_label != "":
                    site3_start = input("\n(FREEDA) (optional) Pick first residue number for region: %s\n" % site3_label)
                    if not re.match(r"^([0-9]+$)", site3_start):
                        print("\nInvalid residue number for region: %s\n" % site3_label)
                    site3_end = input("\n(FREEDA) (optional) Pick last residue number for region: %s\n" % site3_label)
                    if not re.match(r"^([0-9]+$)", site3_end):
                        print("\nInvalid residue number for region: %s\n" % site3_label)

                # update gene parameters
                all_genes_dict[user_gene] = [advanced_option_1, advanced_option_2,
                                                     [site1_label, site1_start, site1_end],
                                                     [site2_label, site2_start, site2_end],
                                                     [site3_label, site3_start, site3_end]]

            # no labels selected needed for this gene
            else:
                print("\nYou did not label any regions in gene: %s\n" % user_gene)
                # gene labels default
                all_genes_dict[user_gene] = [advanced_option_1, advanced_option_2,
                                                 [site1_label, site1_start, site1_end],
                                                 [site2_label, site2_start, site2_end],
                                                 [site3_label, site3_start, site3_end]]


    all_genes = [gene for gene in all_genes_dict if gene != ""]

    # SET BLAST THRESHOLD
    # initial percent identity threshold for blast matches analysis
    if t is None:
        t_initial = 30
    else:
        t_initial = t

    # generate basic folders for input if not present
    folder_generator.generate_basic_folders(wdir)

    # get all species and genome names
    all_genomes = [genome[1] for genome in genomes_preprocessing.get_names(wdir, ref_species,
                                                                           final_excluded_species, subgroup)]

    # ----------------------------------------#
    ######## GET ALL INPUT DATA  ########
    # ----------------------------------------#

    # generate a reference Genome object
    ref_genome_present, ensembl, ref_species, ref_genomes_path, ref_genome_contigs_dict, \
        biotype, all_genes_ensembl, all_gene_ids_ensembl = input_extractor.generate_ref_genome_object(wdir, ref_species)

    # check if provided gene names are present in ensembl object for ref assembly
    if not input_extractor.validate_gene_names(ref_species, all_genes, all_genes_ensembl, ensembl):
        return

    # stop pipeline if the reference genome is absent
    if not ref_genome_present:
        print("\n...FATAL ERROR... : There is no reference genome detected -> exiting the pipeline now...")
        return

    for gene in all_genes:

        print("\n----------- * %s * -----------" % gene)
        # get structure prediction model from AlphaFold
        possible_uniprot_ids = input_extractor.get_uniprot_id(ref_species, gene, ensembl)
        model_seq, uniprot_id = input_extractor.fetch_structure_prediction(wdir, ref_species,
                                                                               gene, possible_uniprot_ids)
        # get sequence input from ensembl
        input_correct, model_matches_input, microexon_present, microexons = input_extractor.extract_input(wdir,
                ref_species, ref_genomes_path, ref_genome_contigs_dict, ensembl, biotype, gene, model_seq, uniprot_id)

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
    print("Checking genome blast databases...")
    tblastn.run_blast(wdir, ref_species, all_genes, final_excluded_species, subgroup)

    # ----------------------------------------#
    ######## RUN EXON FINDING ########
    # ----------------------------------------#

    if exon_extractor.check_blast_output(ref_species, blast_output_path, t_initial, all_genes):
            result_path = exon_extractor.analyze_blast_results(wdir, blast_output_path, ref_species, t_initial,
                        all_genes, all_genomes, final_excluded_species, gui, logging_window, all_genes_dict, subgroup)
    else:
        return

    # ----------------------------------------#
    ######## RUN PAML and PyMOL ########
    # ----------------------------------------#

    # run PAML
    nr_of_species_total_dict, PAML_logfile_name, day, failed_paml, \
    genes_under_pos_sel = paml_launcher.analyze_final_cds(wdir, ref_species, result_path,
                    all_genes, codon_frequencies, gui, logging_window, final_excluded_species, subgroup)

    # visualize PAML result
    paml_visualizer.analyze_PAML_results(wdir, result_path, all_genes,
                            nr_of_species_total_dict, ref_species, PAML_logfile_name, day, genes_under_pos_sel,
                            failed_paml, codon_frequencies, gui, subgroup)
    # run PyMOL
    for gene in all_genes:

        # do not allow further analysis of failed paml runs
        if gene in failed_paml:
            continue

        # check if model seq and input seq match and check if exactly one model exists
        elif structure_builder.check_structure(wdir, ref_species, gene):
            successful = structure_builder.run_pymol(wdir, ref_species, result_path,
                                                         gene, genes_under_pos_sel, all_genes_dict)
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

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-d", "--wdir",
                        help="specify working directory (absolute path to Data folder ex. /Users/user/Data/)", type=str,
                        default=None)
    parser.add_argument("-rs", "--ref_species",
                        help="specify reference organism (e.g. Mm)", type=str, default=None)
    parser.add_argument("-t", "--blast_threshold",
                        help="specify percentage identity threshold for blast (from 30 to 80)", type=int, default=None)
    parser.add_argument("-cf", "--codon_frequencies",
                        help="specify codon frequency models (F3X4 or F61)", type=str, default=None)
    parser.add_argument("-sg", "--subgroup",
                        help="specify subgroup (e.g. hominoidea)", type=str, default=None)
    parser.add_argument("-es", "--excluded_species",
                        help="specify species to exclude (e.g. Ha Gd)", type=str, default=None)

    args = parser.parse_args()
    freeda_pipeline(wdir=args.wdir,
                    ref_species=args.ref_species,
                    t=args.blast_threshold,
                    codon_frequencies=args.codon_frequencies,
                    subgroup=args.subgroup,
                    excluded_species=args.excluded_species)