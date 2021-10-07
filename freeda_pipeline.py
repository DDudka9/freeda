#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:36:16 2021

@author: damian
Main module of the freeda package. Takes user input and performs automatic input extraction, tblastn, exon finding
and molecular evolution analysis (PAML) followed by overlay of putative adaptive sites onto 3D structure (PyMOL).

"""

"""
ISSUE -> Reference genome (mouse) lacks gene name Ap2m1 -> but ensembl has it, model was found so whats the problem? -> gene_name issue?


Traceback (most recent call last):
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda_pipeline.py", line 409, in <module>
    freeda_pipeline(ref_species=args.ref_species, t=args.blast_threshold, wdir=args.wdir)
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda_pipeline.py", line 290, in freeda_pipeline
    ref_species, int(t), all_proteins, all_genomes)
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda/exon_extractor.py", line 59, in analyse_blast_results
    matches = matches_generator.generate_matches(match_path, t, protein_name, genome_name, genome_index)
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda/matches_generator.py", line 39, in generate_matches
    concatenated_matches = split_large_contigs(dataframes).reset_index(drop=True)
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda/matches_generator.py", line 151, in split_large_contigs
    concatenated_matches = pd.concat([i for i in list_new_matches])
  File "/Users/damian/anaconda3/envs/py37/lib/python3.7/site-packages/pandas/core/reshape/concat.py", line 284, in concat
    sort=sort,
  File "/Users/damian/anaconda3/envs/py37/lib/python3.7/site-packages/pandas/core/reshape/concat.py", line 331, in __init__
    raise ValueError("No objects to concatenate")



"""

# TODO:
#    Deal with no matches at threshold (see above)
#    Cpne9 -> Mi -> exon 10 has a small repetition that causes exon splitting (10 and 11) and makes total exon number higher -> CDS not in frame -> all exons are a mess
#    Also exon 22 (the one that should not exist) is called exon 23 ???
#    That is no longer true in Ms !! (exon found perfectly well, except that we lost 3bp from exon 10 that aligned with an intron)
#    I should generally double check exon finding in contigs that show deletions -> might be exactly that same issue
#    SOLUTION : fix the exon finding function -> test on Cpne9 Mi contig QGOO01037269.1__rev
#    BETTER SOLUTION : since this predicts the whole sequence to be out of frame -> get an alignment score for each sequence in final alignment (against ref) -> low score -> remove from alignment -> FIXED?
#    Izumo3 -> Mc -> there is 10bp frameshift insertion in last exon (7) -> FREEDA deletes it in "corrected"
#    -> blastp with uniprot Izumo3 shows that this is consistent with known protein seq and that framshift is most likely a seq error (so good job FREEDA)
#    Ptprd -> 5,6,8 microexons and 500kb gene -> "stich" missing bp in that case as if it was a single microexon
#    Rnf187 -> is misssing START codon (exon 4 -> 3bp; which is a STOP codon) -> that transcript does not have a START codon in ensembl ("START lost")
#    Run PAML on Rnf187 to see how does it affect the pipeline
#    When running Rnf187 I get "data missing" FATAL ERROR -> and there is no "model_matches_input.txt" -> why? -> probably cose no START codon (no cds so no comparison with model)
#    The pipeline will run Rnf187 till the end but it doesnt run from middle (no "model_matches_input.txt" file !!!)
#    Update no matches found at given treshold -> min 7 species need to be present -> otherwise issue WARNING but not FATAL ERROR
#    Refactor : change "protein_name" to "gene_name" and "protein" to "gene_name"
#    Think about PAML visualization and what the "gray" bars mean -> longest sequence, not reference, what about 3D overlay?
#    Use Mo for prediction of sequence accuracy -> compare with NCBI Mo
#    Get full gene name list and pass it to GUI -> user can only pick valid gene names
#    Check if NCBI datasets can give uniprot ID -> is it better than pyensembl?
#    To check operation system -> os.uname().sysname -> macOS is "Darwin", linux is "Linux"
#    Figure out how to bypass the nead for pyensembl install release
#    Issue with excel sheet in Mis18bp1 -> problem with finding coverage? -> FIXED?
#    Issue with PAML visualization graph -> CDS should be that of the reference species but Mis18bp1 it looks like its the longest's species # I decided that its ok
#    Fix the bioservices issue (Brian) -> in virtual box and pyinstaller the colorlog module doesnt have "logging" attribute -> deprecated in python 3.8 ?
#    Use colorlog module to colour the log files
#    TESTING > 10kb flanks on CD46 and CD55 with 70 t and 30kb flanks (08_22_2021)
#                   -> ISSUE -> flanks are as big as the split_large_contigs function -> might be getting same matches on artificially different contigs???
#                               -> requires testing but probably not (CD46)
#           -> CD46 C-terminus was successfully recovered!!! (so larger flanks help -> need to be paired with higher thresholds though)
#           -> Try dynamic flanking -> 10kb if gene < 30kb and 30kb if gene > 30kb
#           -> CD46 ended up NOT passing positive selection tests (LRT 2.24) -> try to run it with species tree? (but the gene tree looks fine)
#           -> try to test flanks 10kb with blastn on CD46 -> NEED TO HAVE CDS IN BLAST INPUT -> it recovers most exons at 30 t but not all (MULATTA 13 exon missing)
#   CONTINUE TESTING -> allowed first exons to be divergent (08_21_2021) -> but it doesnt work -> N-term needs to pass synteny check first
#    0) Use Apbb1 - Ay -> contig LIPJ01008178.1__rev -> exon 11 has one single N and it gets thrown out -> fix conservatively? -> or more conservative would be to delete that base -> gBlocks will take care of the frameshift
#    0) AP2M1 -> cannot overlay on 3D structure cose of microexon but should still show model
#    0) ISSUE with "STOP codon detected in öAST exon (24) in Gorilla Numa1 -> last exon is microexon (25) so its missing but finder thinks there is a STOP in 24 (which there is not)
#           -> also C-term synteny check should not run if last exon is missing (currently exon 24 in Gorilla is syntenic) -> probably DONE
#           -> also add bp number to microexon info in model_incompatible.txt file and log it in exon finder -> DONE
#    1) ISSUE with Haus8 -> Gs -> SRMG01015959.1__for -> part of exon 4 does not align (the other one does), there is insertion as well
#                           it created a frameshift at the beginning of the sequence (22aa) present in translated alignment
#                           it might skew the PAML result for Haus8
#                           SOLUTION : drop Gs from Haus8 analysis
#                           Generally Haus8 rat has some very divergent regions but they match the rat uniprot sequence
#                           Haus8 is a weird protein -> possibly many duplications, retrotranspositions
#    2)  TESTING run time for same protein using higher blast thresholds (50 and 70)
#    3)  ISSUE with CD46 primates -> 10-13 exons found only in mulatta -> check blast file
#                       Consider running a blastn (nucleotide) instead of tblastn (protein) -> tried that, still doesnt find all exons
#                       Consider extending the arms above 10kb to 30kb to check if thats the issue (probably same as CD55)
#    6) ISSUE with defining parameters:
#           Define a module for tweaking parameters (advanced_parameters.py)
#               - duplication restriction (switches on the duplication score)
#               - flanking arms (default 10kb) -> recommend for large introns (ex. primate default to 30kb)
#               - blast threshold
#               - homology threshold
#               - synteny threshold
#               - coverage threshold
#               - non_ACGT corrector (to mirror CDS position)
#               - pymol residues
#               - input known sites
#    7) ISSUE with BEB results for non-adaptive proteins:
#            Something weird about Bub1 -> lots of >0.90 sites but M7 higher than M8
#            Same with Cenp-W
#            Not sure what the solution is -> I made sure proteins that do not score in M8 vs M7 are not visualized
#    8) ISSUE with the cds_cloner function (requires refactoring):
#           Cloner module needs revision to get hamming distance duplication comparison compare
#           the actual duplicated exons and not only the number of exon they carry
#           test on Aurkc Ap
#    9) ISSUE with exon_finding function:
#           Single non_ACGT bases currently lead to whol exon loss
#           SOLUTION: THINK ABOUT FLIPPING non_ACGT INTO CORRESPONDING CDS POSITION (conservative)
#           this could save these exons!
#    10) ISSUE with early SROP codons :
#           THERE IS AN ISSUE WITH: if earlier STOP present in other species then
#           ref species gets translated normally and final_ref_dict is +1
#           which leads to ValueError in get_omegas function
#           SOLUTION: enabled the STOP_remover function( FIXED?)
#           TO FIX: dashes in the MAFFT alignment (need to remove these positions before
#           counting codons -> test on Haus8)
#           08_21_2021 -> I dont really know what this comment mean
#    11) ISSUE with correction:
#           THERE IS AN ERROR IN HAUS8 CORRECTION function -> not same lengths?
#           check the print screen
#    13) ISSUE with the log files:
#           Make FREEDA log file more readable (indentations)
#    14) ISSUE with running Ap2m1:
#           How come "Contig too short to check C-term synteny 0bp aligned" for contig 81143 in genome11 Ap2m1
#           SOLUTION: Probably connected to exon4 being a 6bp microexon and NOT deleted from exon input but
#           There are 13 exons expected instead of 11 -> exon 12 is skipped for some reason; alignment looks good
#           Also alignment of single exons from exon 7 is messed up (linux default file order problem again?)
#           early STOP remover function worked well -> post trimming it was easier to align hence difference in "no_STOP" alignment length
#           but since last 4 single exons were aligned poorly, the stop codons were missing/were displaced in other species
#    16)  TESTING: download genomes of more primates and try reproducing Schuler 2010 MBE paper

print("\nImporting all modules and libraries...\n")

from freeda import input_extractor
from freeda import tblastn
from freeda import exon_extractor
from freeda import paml_launcher
from freeda import paml_visualizer
from freeda import structure_builder
from freeda import genomes_preprocessing
import os


def freeda_pipeline(wdir=None, ref_species=None, t=None):
    # current directory must be the "Data" folder

    if wdir is None:
        wdir = os.getcwd() + "/"

    if ref_species is None:
        ref_species = "Mm"

    # reference species sequences: protein seq, cds, exons, gene (ex. Mus musculus)
    if ref_species != "Mm" and ref_species != "Hs":
        print("\nSupported reference species are mouse : %s and human : %s" % ('"Mm"', '"Hs"'))
        return

    # initial percent identity threshold for blast matches analysis
    if t is None:
        t = 30

    user_input0 = None
    user_input1 = None
    user_input2 = None
    user_input3 = None

    # ----------------------------------------#
    ######## GET USER INPUT ########
    # ----------------------------------------#

    print("Choose which parts of the pipeline you would like to run (all 'y' is a good strategy for single proteins) : ")
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
    input_extractor.generate_basic_folders(wdir)
    # get all proteins to be analysed (these are gene names)
    all_proteins = [protein.rstrip("\n") for protein in open(wdir + "proteins.txt", "r").readlines() if protein != "\n"]


    # get all species and genome names
    all_genomes = [genome[1] for genome in genomes_preprocessing.get_names(ref_species, ref_genome=False)]

    #all_species = [names[0] for names in all_names]
    #all_genome_names = [names[1] for names in all_names]


    # check if the user had previously obtained data for given list of proteins
    if user_input0 == "n":
        input_present = True

        for protein in all_proteins:
            structure_path = wdir + "Structures/" + protein + "_" + ref_species
            if "model_matches_input_seq.txt" in os.listdir(structure_path) or "model_incompatible.txt" in os.listdir(structure_path):
                print("\nAll input data and structure model for : %s are present." % protein)
            else:
                input_present = False
                print("\n...WARNING... : Data are missing for : %s" % protein)
        if not input_present:
            print("\n...FATAL_ERROR... : Input data for at least one protein are missing -> exiting the pipeline now...")
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
        ref_genome_present, ensembl, ref_species, ref_genomes_path, ref_genome_contigs_dict, \
                        biotype, all_genes_ensembl = input_extractor.generate_ref_genome_object(wdir, ref_species)

        # check if provided gene names are present in ensembl object for ref assembly
        if not input_extractor.validate_gene_names(all_proteins, all_genes_ensembl):
            return

        # stop pipeline if the reference genome is absent
        if not ref_genome_present:
            print("\n...FATAL ERROR... : There is no reference genome detected -> exiting the pipeline now...\n"
                  "\n   Make sure you downloaded it into ../Data/Reference_genomes from "
                  " https://www.ncbi.nlm.nih.gov/assembly -> (mouse: GCA_000001635.8; human: GCA_000001405.28) -> "
                  "GenBank -> Genomic FASTA(.fna)")
            return

        for protein in all_proteins:

            print("\n----------- * %s * -----------" % protein)
            # get structure prediction model from AlphaFold
            possible_uniprot_ids = input_extractor.get_uniprot_id(ref_species, protein)
            model_seq, uniprot_id = input_extractor.fetch_structure_prediction(wdir, ref_species, protein, possible_uniprot_ids)
            # get sequence input from ensembl
            input_correct, model_matches_input, microexon_present, microexons = input_extractor.extract_input(
                wdir, ref_species, ref_genomes_path,
                ref_genome_contigs_dict, ensembl, biotype,
                protein, model_seq, uniprot_id
            )

            if input_correct:
                print("\nInput data have been generated for protein: %s\n\n" % protein)

            if not input_correct:
                print("\n...FATAL ERROR... : Input data generation FAILED for protein: %s -> exiting the pipeline now...\n" % protein)
                return

            if not model_matches_input:
                print("...WARNING... : Structure prediction for protein: %s DOES NOT have a match in available ensembl "
                        "database -> cannot overlay FREEDA results onto a 3D structure\n" % protein)
                print("...WARNING... : Protein may still be analyzed using PAML but without 3D structure overlay\n")

            if microexon_present:
                print("...WARNING... : Sequence for: %s found in Ensembl contains microexons : %s\n" % (protein, microexons))
                print("...WARNING... : Microexons are difficult to align and are removed\n")

    # ----------------------------------------#
    ######## RUN BLAST ########
    # ----------------------------------------#

    blast_output_path = wdir + "Blast_output/"
    if user_input1 == "y":
        print("Checking genome blast databases...")
        tblastn.run_blast(wdir, ref_species, all_proteins)
        #if blast_output_path is None:
        #    print("\n...FATAL ERROR... : Blast database build failed for at least one genome"
        #          "\n                               -> exiting the pipeline now...")
        #return

    # ----------------------------------------#
    ######## RUN EXON FINDING ########
    # ----------------------------------------#

    if user_input2 == "y":
        if exon_extractor.check_blast_output(blast_output_path, t, all_proteins):
            result_path = exon_extractor.analyse_blast_results(wdir, blast_output_path,
                                                               ref_species, int(t), all_proteins, all_genomes)
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
                nr_of_species_total_dict, \
                PAML_logfile_name, day, \
                failed_paml, \
                proteins_under_positive_selection = paml_launcher.analyse_final_cds(wdir, ref_species,
                                                                                    result_path, all_proteins)

                if not all([nr_of_species_total_dict]):
                    print("\n...FATAL_ERROR... : Failed PAML analysis -> exiting the pipeline now ...")
                    return

                # visualize PAML result
                paml_visualizer.analyse_PAML_results(wdir, result_path, all_proteins,
                                                     nr_of_species_total_dict, ref_species,
                                                     PAML_logfile_name, day, proteins_under_positive_selection)
                # run PyMOL
                for protein in all_proteins:
                    # do not allow further analysis of failed paml runs

                    if protein in failed_paml:
                        continue
                    # check if model seq and input seq match and check if exactly one model exists
                    elif structure_builder.check_structure(wdir, ref_species, protein):
                        successful = structure_builder.run_pymol(wdir, ref_species, result_path,
                                                                 protein, proteins_under_positive_selection,
                                                                 offset=None)
                        if not successful:
                            print("\nThe structure for : %s was not built successfully." % protein)
                            continue
                    else:
                        print("\nPrediction model for : %s DOES NOT match input sequence "
                              "-> cannot overlay FREEDA results onto a 3D structure\n" % protein)

    if user_input3 == "y" and user_input2 == "y":

        # run PAML
        nr_of_species_total_dict, \
        PAML_logfile_name, day, \
        failed_paml, \
        proteins_under_positive_selection = paml_launcher.analyse_final_cds(wdir, ref_species, result_path, all_proteins)

        if not all([nr_of_species_total_dict]):
            print("\n...FATAL_ERROR... : Failed PAML analysis -> exiting the pipeline now ...")
            return

        # visualize PAML result
        paml_visualizer.analyse_PAML_results(wdir, result_path, all_proteins,
                                             nr_of_species_total_dict, ref_species,
                                             PAML_logfile_name, day, proteins_under_positive_selection)
        # run PyMOL
        for protein in all_proteins:
            # do not allow further analysis of failed paml runs

            if protein in failed_paml:
                continue
            # check if model seq and input seq match and check if exactly one model exists
            elif structure_builder.check_structure(wdir, ref_species, protein):
                successful = structure_builder.run_pymol(wdir, ref_species, result_path,
                                                         protein, proteins_under_positive_selection,
                                                         offset=None)
                if not successful:
                    print("\nThe structure for : %s was not built successfully." % protein)
                    continue
            else:
                print("\nPrediction model for : %s DOES NOT match input sequence "
                      "-> cannot overlay FREEDA results onto a 3D structure\n" % protein)

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
                        help="specify reference organism (default is mouse)", type=str, default="Mm")
    parser.add_argument("-t", "--blast_threshold",
                        help="specify percentage identity threshold for blast (default is 30)", type=int, default=70)


    args = parser.parse_args()
    freeda_pipeline(ref_species=args.ref_species, t=args.blast_threshold, wdir=args.wdir)



"""

 
                # get structure model using AlphaFold url request
                prediction_url = input_extractor.get_prediction(wdir, ref_species, protein)
                if prediction_url == None:
                    print("AlphaFold prediction not available for: %s\n" % protein)
                    model_equal_input = False
                    pass
                elif prediction_url == True:
                    print("Structure prediction model for: %s already exists\n" % protein)
                    # check if model sequence equals the blasted sequence
                    model_equal_input = structure_builder.compare_model_with_input(wdir, ref_species, protein)
                    pass
                else:
                    print("\n(FREEDA) Please input structure prediction for protein: %s\n(copy the following url into your browser " \
                      "-> click PDB file -> save in ../Data/Structures/%s)\n\n " \
                         "%s\n\n ...WARNING... Verify protein identity (if incorrect find model in AlphaFold browser)" 
                                                     % (protein, protein + "_" + ref_species, prediction_url))
                    input("\n(FREEDA) When done press ENTER\n")
                    # check if model sequence equals the blasted sequence
                    model_equal_input = structure_builder.compare_model_with_input(wdir, ref_species, protein)


"""
