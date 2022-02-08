#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:36:16 2021

@author: damian
Main module of the freeda package. Takes user input and performs automatic input extraction, tblastn, exon finding
and molecular evolution analysis (PAML) followed by overlay of putative adaptive sites onto 3D structure (PyMOL).

"""

"""

Traceback (most recent call last):
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda_pipeline.py", line 378, in <module>
    freeda_pipeline(ref_species=args.ref_species, t=args.blast_threshold, wdir=args.wdir)
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda_pipeline.py", line 310, in freeda_pipeline
    gene, genes_under_pos_sel)
TypeError: run_pymol() missing 1 required positional argument: 'all_genes_dict'

Process finished with exit code 1


"""

# TODO:
#    0) ESSENTIAL -> test CENPP in Cf and check compatibility cose Fc doesnt have it annotated
#    0) ESSENTIAL -> cloned cds frameshift check should not penalize gaps or at least not say "frameshift"
#                               because Cj in NRLP11 is called frameshift despite having just an exon missing
#                                    -> I have to rerun it with lower threshold because Cj did not blast enough
#    15) ISSUE with early STOP codons :
#           THERE IS AN ISSUE WITH: if earlier STOP present in other species then
#           ref species gets translated normally and final_ref_dict is +1
#           which leads to ValueError in get_omegas function
#           SOLUTION: enabled the STOP_remover function( FIXED?)
#           TO FIX: dashes in the MAFFT alignment (need to remove these positions before
#           counting codons -> test on Haus8)
#           08_21_2021 -> I dont really know what this comment mean

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


def freeda_pipeline(wdir=None, ref_species=None, t=None, codon_frequencies=None):
    """Main function running all freeda pipeline from command line"""

    if wdir is None:
        wdir = os.getcwd() + "/"

    if ref_species is None:
        ref_species = "Mm"

    # initial percent identity threshold for blast matches analysis
    if t is None:
        t = 30

    if codon_frequencies is None:
        codon_frequencies = "F3X4"

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

    # get settings
    aligner = "mafft"

    # get all species and genome names
    all_genomes = [genome[1] for genome in genomes_preprocessing.get_names(wdir, ref_species, ref_genome=False)]

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
    ######## GET ALL INPUT DATA  ########
    # ----------------------------------------#

    # User wants to generate input data
    if user_input0 == "y":

        # generate a reference Genome object
        # WHY DOES IT RETURN REF_SPECIES?
        ref_genome_present, ensembl, ref_species, ref_genomes_path, ref_genome_contigs_dict, \
                        biotype, all_genes_ensembl = input_extractor.generate_ref_genome_object(wdir, ref_species)

        # check if provided gene names are present in ensembl object for ref assembly
        if not input_extractor.validate_gene_names(all_genes, all_genes_ensembl):
            return

        # stop pipeline if the reference genome is absent
        if not ref_genome_present:
            print("\n...FATAL ERROR... : There is no reference genome detected -> exiting the pipeline now...\n"
                  "\n   Make sure you downloaded it into ../Data/Reference_genomes from "
                  " https://www.ncbi.nlm.nih.gov/assembly -> (mouse: GCA_000001635.8; human: GCA_000001405.28) -> "
                  "GenBank -> Genomic FASTA(.fna)")
            return

        # get names of

        for gene in all_genes:

            print("\n----------- * %s * -----------" % gene)
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
                print("\nInput data have been generated for gene: %s\n\n" % gene)

            if not input_correct:
                print("\n...FATAL ERROR... : Input data generation FAILED for gene: %s - please remove from analysis"
                      " -> exiting the pipeline now...\n" % gene)
                return

            if not model_matches_input:
                print("...WARNING... : No matching structure prediction model is available for : %s "
                      "-> cannot overlay FREEDA results onto a 3D structure\n" % gene)
                print("...WARNING... : gene will still be analyzed using PAML but without 3D structure overlay\n")


    # ----------------------------------------#
    ######## RUN BLAST ########
    # ----------------------------------------#

    blast_output_path = wdir + "Blast_output/"
    if user_input1 == "y":
        print("Checking genome blast databases...")
        tblastn.run_blast(wdir, ref_species, all_genes)

    # ----------------------------------------#
    ######## RUN EXON FINDING ########
    # ----------------------------------------#

    if user_input2 == "y":
        if exon_extractor.check_blast_output(blast_output_path, t, all_genes):
            result_path = exon_extractor.analyze_blast_results(wdir, blast_output_path,
                                                    ref_species, int(t), all_genes, all_genomes, aligner)
        else:
            print("\n     ...FATAL ERROR... : Genome of at least one species contains "
                  "no matches above the identity threshold used : %s -> use a lower one " 
                    "-> exiting the pipeline now..." % t)
            return

    # ----------------------------------------#
    ######## RUN PAML and PyMOL ########
    # ----------------------------------------#

    if user_input3 == "y" and user_input2 == "n":

        nr_of_tries = 1
        while nr_of_tries <= 3:
            user_input4 = input("----->  Indicate folder with results (no slashes or quotes): ")
            result_path = wdir + user_input4 + "/"

            if os.path.isdir(result_path) is False:
                print("\n(FREEDA) I could not find your results folder (%s/3)" % nr_of_tries)
                nr_of_tries += 1
            else:
                nr_of_tries = float("inf")
                result_path = wdir + user_input4 + "/"

                # run PAML
                nr_of_species_total_dict, PAML_logfile_name, day, \
                failed_paml, genes_under_pos_sel = paml_launcher.analyze_final_cds(wdir, ref_species,
                                                            result_path, all_genes, aligner, codon_frequencies)

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
                                                            result_path, all_genes, aligner, codon_frequencies)

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
                        help="specify percentage identity threshold for blast (default is 60)", type=int, default=60)
    parser.add_argument("-f", "--codon_frequencies",
                        help="specify codon frequency models (F3x4 is default)", type=str, default="F3x4")

    args = parser.parse_args()
    freeda_pipeline(wdir=args.wdir, ref_species=args.ref_species, t=args.blast_threshold,
                    codon_frequencies=args.codon_frequencies)

