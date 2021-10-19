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

SUGGESTIONS:
#               Use Mo for prediction of sequence accuracy -> compare with NCBI Mo
#               Bioservices 1.8.1 release is ready -> not sure if I need to update from 1.7.12
#               Use pytest -> make TestClass for all tests and use "test_*.py" notation for the testing module
#               Check if NCBI datasets can give uniprot ID -> is it better than pyensembl?
#               To check operation system -> os.uname().sysname -> macOS is "Darwin", linux is "Linux"
#               Figure out how to bypass the nead for pyensembl install release
#               Use colorlog module to colour the log files
#               Get full gene name list and pass it to GUI -> user can only pick valid gene names
#               CONTINUE TESTING -> allowed first exons to be divergent (08_21_2021) -> but it doesnt work -> N-term needs to pass synteny check first


----------- * Nudt11 * -----------

Looking for structure in AlphaFold database for: Nudt11 ...

Based on uniprot id: P0C028 structure prediction for protein: Nudt11 has been found (164 aa)

Exons FAILED to assemble expected CDS for: Nudt11-201


...FATAL ERROR... : Input data generation FAILED for protein: Nudt11 -> exiting the pipeline now...


-> its a protein coding gene with 2 exons where exon 2 is a SINGLE nucleaotide (A) - part of the STOP









 --------- * Prdm9 * --------- 


Alignment score for species : >Mm = -0.8530805687203792
...WARNING... : CDS for species : >Mm in Prdm9 protein contains a frameshift -> eliminated from alignment

Alignment score for species : >Ay = -1.0793838862559242
...WARNING... : CDS for species : >Ay in Prdm9 protein contains a frameshift -> eliminated from alignment

Alignment score for species : >An = -1.0209320695102686
...WARNING... : CDS for species : >An in Prdm9 protein contains a frameshift -> eliminated from alignment

Alignment score for species : >Rs = 0.872827804107425
Alignment score for species : >Rn = -1.032780410742496
...WARNING... : CDS for species : >Rn in Prdm9 protein contains a frameshift -> eliminated from alignment

 
 WARNING : Insertions in positions in ref MSA deleted: [1491, 1492, 1493, 1958, 3698, 3699]
 
---- Insertions detected in Prdm9 alignment -> these bp positions were removed in all species forcing conserved alignment

[]

...WARNING... : Failed PAML analysis for : Prdm9 -> probably not enough species in the alignment (e.g. poor alignment of repetitive regions)

 --------------->  PAML analysis completed in 0.01356441577275594 minutes or 0.00022607386112213134 hours
Traceback (most recent call last):
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda_pipeline.py", line 434, in <module>
    freeda_pipeline(ref_species=args.ref_species, t=args.blast_threshold, wdir=args.wdir)
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda_pipeline.py", line 356, in freeda_pipeline
    day, prots_under_pos_sel, failed_paml)
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda/paml_visualizer.py", line 57, in analyse_PAML_results
    get_alignment_matching_structure(result_path, ref_species, protein, all_matched_adaptive_sites_ref)
  File "/Users/damian/PycharmProjects/freeda_2.0/freeda/paml_visualizer.py", line 75, in get_alignment_matching_structure
    aa_features_dict = dictionary[protein]
KeyError: 'Prdm9'


"""

# TODO:
#    0) ESSENTIAL -> Should all 17 + 1 rodent genomes be used? Some may introduce problematic alignments
#                           and mixed duplications events -> Ha, Pd, Mn, Gd, Ap, Ay
#    0) ESSENTIAL -> what to call adaptively evolving? -> both M1a vs M2a and M7 vs M8 should be < 0.05 ?
#                   -> e.g. Cxxc1 is unlinkely rapidly evolving but it scores in M7 vs M8
#    0) ESSENTIAL -> Sgo2b last exon (7) is called "None" in exon finder in all contigs -> last bp in cds is missing in gene
#           SOLUTION : Special case for last bp?
#    0) ESSENTIAL -> Sgo2b has a frameshift deletion in exon 6 -> freeda makes it inf and takes Sgo2a as true Sgo2b
#           SOLUTION : Deactivate frameshift check? Sometimes frameshifts might be real
#    0) ESSENTIAL -> Prdm9 fails cose of few sequences passing 90% -> cose last exon (10) is huge and 8kb distant -> and there is up to 5000 matches from blast
#           so probably the final threshold of 65 misses some matches in exon 10 and contig ends early -> also plenty of Ns in exon 10 in Mi, Mc etc.
#           Also exon 10 in Prdm9 has duplications that make alignment difficult -> thats why matches and locus not always align (probably main culprit)
#           Main problem -> Mafft seems to fail aligning Prdm9_final.fasta (but Ugene mafft does it with no problem)
#           SOLUTIONS : 1) Nevertheless FREEDA should not crash; 2) Increase matches limit from 40 to 60 and run again; 3) Decrease threshold from 90% to 70% and run again
#           I might need to add single exon length checkpoint -> dash the exon that is double the size of the ref exon -> DONE
#           SOLVED -> exon insertion checkpoint at single exon mapping saved Prdm9 aligmnent and showed rapid evolution of DNA binding repeats in exon 10
#    0) UPGRADE -> Visualization module should get info about domains as well
#    0) UPGRADE -> Interpro module often overlaps domains (ex. Dmc1 -> 4 overlapping domains)
#                       -> if coordinates overlap > 50% -> pick only first domain -> FIXED?
#    0) ESSENTIAL -> test using different aligners - not for user - (Clustal Omega, Muscle, PRANK)
#    0) ESSENTIAL -> run proteins that have seqs from Mo on uniprot -> compare
#    1) ISSUE  -> Ask Tim from pyensembl how to get release outside command line -> DONE?
#    2) ISSUE  -> Ptprd -> 5,6,8 microexons and 500kb gene -> "stich" missing bp in that case as if it was
#                  a single microexon
#                  SOLUTION : no good solution for that so far, stiching exons does not help cose of flanking exons
#    3) ISSUE   -> Rnf187 -> is misssing START codon (exon 4 -> 3bp; which is a STOP codon)
#                   -> that transcript does not have a START codon in ensembl ("START lost")
#                   Run PAML on Rnf187 to see how does it affect the pipeline
#                   When running Rnf187 I get "data missing" FATAL ERROR -> and there is no "model_matches_input.txt"
#                   -> why? -> probably cose no START codon (no cds so no comparison with model)
#                   The pipeline will run Rnf187 till the end but it doesnt run from middle
#                   (no "model_matches_input.txt" file !!!)
#    3) ISSUE   -> Refactor : change "protein_name" to "gene_name" and "protein" to "gene_name"
#    4) ISSUE   -> Fix the bioservices issue (Brian) -> in virtual box and pyinstaller the colorlog module
#                   doesnt have "logging" attribute -> deprecated in python 3.8 ?
#    5) ISSUE   -> flanks are as big as the split_large_contigs function -> might be getting same matches
#                      on artificially different contigs??? -> requires testing but probably not (CD46)
#               -> CD46 C-terminus was successfully recovered!!! (so larger flanks help -> need to be
#                           paired with higher thresholds though)
#               -> Try dynamic flanking -> 10kb if gene < 30kb and 30kb if gene > 30kb
#               -> CD46 ended up NOT passing positive selection tests (LRT 2.24) -> try to run it with species tree?
#                               (but the gene tree looks fine)
#               -> try to test flanks 10kb with blastn on CD46 -> NEED TO HAVE CDS IN BLAST INPUT
#                           -> it recovers most exons at 30 t but not all (MULATTA 13 exon missing)
#    6) Use Apbb1 - Ay -> contig LIPJ01008178.1__rev -> exon 11 has one single N and it gets thrown out
#                           -> fix conservatively? -> or more conservative would be to delete that base
#                           -> gBlocks will take care of the frameshift
#    7) AP2M1 -> cannot overlay on 3D structure cose of microexon but should still show model
#    8) ISSUE with "STOP codon detected in LAST exon (24) in Gorilla Numa1 -> last exon is microexon (25)
#                   so its missing but finder thinks there is a STOP in 24 (which there is not)
#           -> also C-term synteny check should not run if last exon is missing (currently exon 24
#                           in Gorilla is syntenic) -> probably DONE
#           -> also add bp number to microexon info in model_incompatible.txt file and log it in exon finder -> DONE
#    9)  TESTING run time for same protein using higher blast thresholds (50 and 70)
#    10)  ISSUE with CD46 primates -> 10-13 exons found only in mulatta -> check blast file
#                       Consider running a blastn (nucleotide) instead of tblastn (protein)
#                                -> tried that, still doesnt find all exons
#                       Consider extending the arms above 10kb to 30kb to check if thats the issue
#                               (probably same as CD55)
#    11) ISSUE with defining parameters:
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
#    12) ISSUE with BEB results for non-adaptive proteins:
#            Something weird about Bub1 -> lots of >0.90 sites but M7 higher than M8
#            Same with Cenp-W
#            Not sure what the solution is -> I made sure proteins that do not score in M8 vs M7 are not visualized
#    13) ISSUE with the cds_cloner function (requires refactoring):
#           Cloner module needs revision to get hamming distance duplication comparison compare
#           the actual duplicated exons and not only the number of exon they carry
#           test on Aurkc Ap
#           test on Nlrp5 Ap BDUI01009531.1__rev
#    14) ISSUE with exon_finding function:
#           Single non_ACGT bases currently lead to whole exon loss
#           SOLUTION: THINK ABOUT FLIPPING non_ACGT INTO CORRESPONDING CDS POSITION (conservative)
#           this could save these exons! -> better woould be to just delete that base and let frameshit be taken care of by Gblocks
#           this function should target only single base non_ACGT and flag these exons -> test on Cxxc1 Ap aligned_rev_comp_BDUI01029349.1__rev exon 7
#    15) ISSUE with early STOP codons :
#           THERE IS AN ISSUE WITH: if earlier STOP present in other species then
#           ref species gets translated normally and final_ref_dict is +1
#           which leads to ValueError in get_omegas function
#           SOLUTION: enabled the STOP_remover function( FIXED?)
#           TO FIX: dashes in the MAFFT alignment (need to remove these positions before
#           counting codons -> test on Haus8)
#           08_21_2021 -> I dont really know what this comment mean
#    16) ISSUE with correction:
#           THERE IS AN ERROR IN HAUS8 CORRECTION function -> not same lengths?
#           check the print screen
#    18) ISSUE with running Ap2m1:
#           How come "Contig too short to check C-term synteny 0bp aligned" for contig 81143 in genome11 Ap2m1
#           SOLUTION: Probably connected to exon4 being a 6bp microexon and NOT deleted from exon input but
#           There are 13 exons expected instead of 11 -> exon 12 is skipped for some reason; alignment looks good
#           Also alignment of single exons from exon 7 is messed up (linux default file order problem again?)
#           early STOP remover function worked well -> post trimming it was easier to align hence difference in "no_STOP" alignment length
#           but since last 4 single exons were aligned poorly, the stop codons were missing/were displaced in other species

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


    # get settings
    aligner = "mafft"


    # get all species and genome names
    all_genomes = [genome[1] for genome in genomes_preprocessing.get_names(ref_species, ref_genome=False)]

    #all_species = [names[0] for names in all_names]
    #all_genome_names = [names[1] for names in all_names]


    # check if the user had previously obtained data for given list of proteins
    if user_input0 == "n":
        input_present = True

        for protein in all_proteins:
            structure_path = wdir + "Structures/" + protein + "_" + ref_species
            if "model_matches_input_seq.txt" in os.listdir(structure_path) \
                    or "model_incompatible.txt" in os.listdir(structure_path):
                print("\nAll input data and structure model for : %s are present." % protein)
            else:
                input_present = False
                print("\n...WARNING... : Data are missing for : %s" % protein)
        if not input_present:
            print("\n...FATAL_ERROR... : Input data for at least one protein are missing "
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
            model_seq, uniprot_id = input_extractor.fetch_structure_prediction(wdir, ref_species,
                                                                               protein, possible_uniprot_ids)
            # get sequence input from ensembl
            input_correct, model_matches_input, microexon_present, microexons = input_extractor.extract_input(
                wdir, ref_species, ref_genomes_path,
                ref_genome_contigs_dict, ensembl, biotype,
                protein, model_seq, uniprot_id
            )

            if input_correct:
                print("\nInput data have been generated for protein: %s\n\n" % protein)

            if not input_correct:
                print("\n...FATAL ERROR... : Input data generation FAILED for protein: %s "
                      "-> exiting the pipeline now...\n" % protein)
                return

            if not model_matches_input:
                print("...WARNING... : Structure prediction for protein: %s DOES NOT have a match in available ensembl "
                        "database -> cannot overlay FREEDA results onto a 3D structure\n" % protein)
                print("...WARNING... : Protein may still be analyzed using PAML but without 3D structure overlay\n")

            #if microexon_present:
            #    print("...WARNING... : Sequence for: %s found in Ensembl contains microexons : %s\n"
            #          % (protein, microexons))
            #    print("...WARNING... : Microexons are difficult to align and are removed\n")

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
                                                    ref_species, int(t), all_proteins, all_genomes, aligner)
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
                failed_paml, prots_under_pos_sel = paml_launcher.analyse_final_cds(wdir, ref_species,
                                                                                  result_path, all_proteins, aligner)

                #if not all([nr_of_species_total_dict]):
                #    print("\n...FATAL_ERROR... : Failed PAML analysis -> exiting the pipeline now ...")
                #    return

                # visualize PAML result
                paml_visualizer.analyse_PAML_results(wdir, result_path, all_proteins,
                                                     nr_of_species_total_dict, ref_species, PAML_logfile_name,
                                                     day, prots_under_pos_sel, failed_paml)
                # run PyMOL
                for protein in all_proteins:
                    # do not allow further analysis of failed paml runs

                    if protein in failed_paml:
                        continue
                    # check if model seq and input seq match and check if exactly one model exists
                    elif structure_builder.check_structure(wdir, ref_species, protein):
                        successful = structure_builder.run_pymol(wdir, ref_species, result_path,
                                                                 protein, prots_under_pos_sel,
                                                                 offset=None)
                        if not successful:
                            print("\nThe structure for : %s was not built successfully." % protein)
                            continue
                    else:
                        print("\nPrediction model for : %s DOES NOT match input sequence "
                              "-> cannot overlay FREEDA results onto a 3D structure\n" % protein)

    if user_input3 == "y" and user_input2 == "y":

        # run PAML
        nr_of_species_total_dict, PAML_logfile_name, day, \
        failed_paml, prots_under_pos_sel = paml_launcher.analyse_final_cds(wdir, ref_species,
                                                                          result_path, all_proteins, aligner)

        #if not all([nr_of_species_total_dict]):
        #    print("\n...FATAL_ERROR... : Failed PAML analysis -> exiting the pipeline now ...")
        #    return

        # visualize PAML result
        paml_visualizer.analyse_PAML_results(wdir, result_path, all_proteins,
                                             nr_of_species_total_dict, ref_species, PAML_logfile_name,
                                             day, prots_under_pos_sel, failed_paml)
        # run PyMOL
        for protein in all_proteins:
            # do not allow further analysis of failed paml runs

            if protein in failed_paml:
                continue
            # check if model seq and input seq match and check if exactly one model exists
            elif structure_builder.check_structure(wdir, ref_species, protein):
                successful = structure_builder.run_pymol(wdir, ref_species, result_path,
                                                         protein, prots_under_pos_sel,
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
                        help="specify percentage identity threshold for blast (default is 30)", type=int, default=50)

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
