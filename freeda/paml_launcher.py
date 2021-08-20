#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 10:51:27 2021

@author: damian
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 21:48:45 2021

@author: damian

Analyses the final cds, gets a gene tree based on translated cds and runs PAML.

"""

from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import datetime
import time
import glob
import os
import logging
import shutil
import subprocess

from freeda import fasta_reader


def analyse_final_cds(wdir, original_species, result_path, all_proteins):
    
    start_time = time.time()
    day = datetime.datetime.now().strftime("-%m-%d-%Y-%H-%M")
        
    # initiate log file to record PAML analysis by reseting the handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    PAML_logfile_name = "PAML" + day + ".log"
    logging.basicConfig(filename=PAML_logfile_name, level=logging.INFO, format="%(message)s")
    
    # make empty dataframe to store PAML results
    #PAML_df = pd.DataFrame({"Protein Name" : \
    #                        "p-value" : \
    #                        "M2a vs M1a"})
    
    # make an empty template for all possible species
    import copy
    all_species = {}
    with open("species.txt", "r") as f:
        all_species[original_species] = ""
        for species in f.readlines():
            all_species[species.rstrip("\n")] = ""
    
    nr_of_species_total_dict = {} 

    for protein in all_proteins:
        
        # check it this protein was already analysed
        if os.path.isdir(result_path + protein + "/" + "PAML_" + protein):
            
           message = "\n################\n\n PAML analysis has been already performed for : %s (skipping)" % protein
           print(message)
           logging.info(message)
           
           # get how many species were analysed
           path_to_final_species = result_path + protein + "/" + protein + "_final.fasta"
           with open(path_to_final_species) as f:
               nr_of_species = f.read().count(">")
               nr_of_species_total_dict[protein] = nr_of_species
           
        # otherwise proceed with the analysis
        else:
        
            # need to use deepcopy function to make an actual dictionary copy
            final_species = copy.deepcopy(all_species)
            final_species_headers = []
        
            # read each cds fasta file and put them into the empty final_species dict
            cds_file = protein + ".fasta"
        
            if os.path.isfile(cds_file) == False:
                message = "\nThis file has not been found : %s." % cds_file
                print(message)
                logging.info(message)
                continue
        
            with open(cds_file, "r") as f:
                file = f.readlines()
            
                # count lines
                line_nr = 0
                for line in file:
                    line_nr += 1
                
                    # first line is always the original species header
                    if line_nr == 1:
                        original_head = line.split(protein + "_")[1].rstrip("\n")
                        continue
                
                    # second line is always the original species sequence
                    if line_nr == 2:
                        # remove STOP codon if present
                        # original_seq = STOP_remover(line.rstrip("\n"))
                        original_seq = line.rstrip("\n")
                        final_species[original_head] = original_seq
                        continue
                
                    # next lines are headers and sequences of cloned cds
                    if line.startswith(">"):
                        head = line.split(protein + "_")[1].rstrip("\n")
                        continue
                
                    else:
                        # remove STOP codon if present
                        # seq = STOP_remover(line.rstrip("\n"))
                        seq = line.rstrip("\n")
                        seq_no_dashes = seq.replace("-", "")
                        if len(seq_no_dashes) / len(original_seq) >= 0.90: 
                            final_species[head] = seq
                            final_species_headers.append(head)
               
            # log species
            final_species_headers.insert(0, original_species)
            nr_of_species = len((final_species_headers))
            nr_of_species_total_dict[protein] = nr_of_species
            message = "\n --------- * %s * --------- \n\n Final species cloned and aligned (+ original) for %s : %s %s" \
                % (protein, protein, str(nr_of_species), str(final_species_headers))
            print(message)
            logging.info(message)
        
            # write all sequences into a file by reading through the dictionary
            final_cds_file = protein + "_final.fasta"
            with open(final_cds_file, "w") as f:
                for species in final_species:
                    
                    if final_species[species] == "":
                        continue
                    else:
                        f.write(">" + species + "\n")
                        f.write(final_species[species] + "\n")
        
        
            protein_folder_path = result_path + protein
            
            shutil.move(protein + ".fasta", protein_folder_path)
            shutil.move(final_cds_file, protein_folder_path)
        
            # generate a PAML folder for a given protein
            PAML_path = protein_folder_path + "/PAML_" + protein
            os.makedirs(PAML_path)
        
            # copy control_file from working directory into protein_folder_path
            control_file = "control_file.ctl"
            shutil.copy("control_file.ctl", PAML_path + "/" + control_file)
        
            # align the final cds sequences
            out_MAFFT = align_final_cds(protein, final_cds_file, result_path) 
            shutil.move(out_MAFFT, protein_folder_path)

            # check and eliminate insertions that cause dashes in original species 
            correction, corrected_filename = eliminate_all_insertions(protein_folder_path, out_MAFFT)
            if correction == True:
                side_note = " \n---- Insertions detected in " + protein + " alignment -> " \
                    "these bp positions were removed in all species forcing conserved alignment\n"
                print(side_note)
                logging.info(side_note)
                no_dashes_out_MAFFT = corrected_filename
            else:
                no_dashes_out_MAFFT = out_MAFFT
        
            # remove STOP codons
            final_cds_file_no_STOP = STOP_remover(protein_folder_path, no_dashes_out_MAFFT, protein)
            shutil.move(final_cds_file_no_STOP, protein_folder_path)
        
            # run gBLOCK
            out_Gblocks = run_Gblocks(final_cds_file_no_STOP, protein, result_path)
            shutil.move(out_Gblocks, protein_folder_path)
            
            # double check for artificial STOP codons introduced by Gblocks -> force conserved alignment
            # it overwrites the previous out_Gblocks file (same name)
            post_Gblocks_STOP_remover(wdir, protein_folder_path, out_Gblocks, protein)
            
            # translate all sequences (inside of the protein folder)
            translated_path = translate_Gblocks(wdir, protein_folder_path, out_Gblocks, protein, original_species)
        
            # run seqret (file is already in protein folder)
            phylip_path = run_seqret(protein, protein_folder_path, out_Gblocks)
            shutil.copy(phylip_path, PAML_path + "/input.phy")
            
            # run RAxML (and move all the RAxML files to protein folder)
            best_tree_path = run_RAxML(protein, protein_folder_path, translated_path)
            shutil.copy(best_tree_path, PAML_path + "/gene.tree")
        
            # run PAML
            message = "\n.........Running PAML for protein: %s.........\n" % protein
            print(message)
            
            M2a_M1a, M8_M7 = run_PAML(wdir, protein, PAML_path, control_file)
            message = "\n -> PAML p-values for protein %s : M2a v M1a - %s and M8 v M7 - %s" \
                % (protein, str(M2a_M1a), str(M8_M7))
            print(message)
            logging.info(message)
    
    # THERE IS NO POINT OF MARKING THIS DUPLICATION???? EACH TIME ITS DUPLICATED
    
    #if os.path.exists(PAML_logfile_name) == True:
    #    message = "\n(FREEDA) This log file: %s already exists." % PAML_logfile_name
    #    print(message)
    #    logging.info(message)
    #    new_PAML_logfile_name = PAML_logfile_name.rstrip(".log") + "_duplicated.log"
    #    os.rename(PAML_logfile_name, new_PAML_logfile_name)
    shutil.move(wdir + PAML_logfile_name, result_path)

    #else:
    #    shutil.move(wdir + PAML_logfile_name, result_path)
    
    
    # for now final analysis complete message logs into the PAML log file -> change that later
    
    # mark the end of the analysis
    message = ("\n --------------->  PAML analysis completed in %s minutes or %s hours" % \
               ((time.time() - start_time)/60,
                 (time.time() - start_time)/60/60))
    print(message)
    logging.info(message)
    
    return nr_of_species_total_dict, PAML_logfile_name, day


def eliminate_all_insertions(protein_folder_path, out_MAFFT):
    import re
    
    correction = False
    
    alignment = AlignIO.read(protein_folder_path + "/" + out_MAFFT, "fasta")
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]
    
    # check if original species contains dashes
    if "-" in seqs[0][1]:
        
        positions = [p.start() for p in re.finditer("-", seqs[0][1])]
        
        # if yes, make a new file for a corrected alignment
        correction = True
        corrected_filename = "corrected_" + out_MAFFT
        file_path = protein_folder_path + "/" + corrected_filename
        with open(file_path, "w") as f:
        
            # loop through each sequence in the alignment
            for spec, seq in seqs:
                # generate a dict with positions as keys and bp as values (starts at 1)
                d = {}
                for i in range(len(seq)):
                    d[i+1] = seq[i]
                
                # eliminate these positions from given seq
                for position in positions:
                    d[position+1] = ""
                
                # write the corrected seqeunce into a new, corrected file
                new_seq = "".join([a for a in d.values()])
                f.write(">" + spec[2:] + "\n")
                f.write(new_seq + "\n")
            
            # these positions are confusingly printing index+1 (change that -> make positions start with 1)
            
        side_note = " \n WARNING : Insertions in positions in original MSA deleted: %s" % str(positions)
        print(side_note)
        logging.info(side_note)
                    
    # if there are no dashes in original species sequence
    else:
        return correction, out_MAFFT
    
    return correction, corrected_filename


# MODIFY THE INSERTION FUNCTION TO ELIMINATE ALSO 3% != 0 insertions > 1 (mostly artificial insertions)



def run_PAML(wdir, protein, PAML_path, control_file):

    from Bio.Phylo.PAML import codeml
    from Bio.Phylo.PAML.chi2 import cdf_chi2
    
    # change working directory to given PAML protein folder
    os.chdir(PAML_path)
    
    # run PAML
    cml = codeml.Codeml(alignment="input.phy", tree="gene.tree", \
                        out_file="output_PAML", working_dir=PAML_path) 
    cml.run(control_file)
    results = codeml.read("output_PAML")
    ns_sites = results.get("NSsites")
    
    M1a = ns_sites.get(1)
    M1a_lnL = M1a.get("lnL")
    
    M2a = ns_sites.get(2)
    M2a_lnL = M2a.get("lnL")
    
    M7 = ns_sites.get(7)
    M7_lnL = M7.get("lnL")
    
    M8 = ns_sites.get(8)
    M8_lnL = M8.get("lnL")
    
    # there are 2 free parameters in both comparisons
    df_M8_M7 = 2
    df_M2a_M1a = 2
    
    # determine the LRTs (Likelihood Ratio Tests) -> needs to be a positive value
    LRT1 = 2*(M2a_lnL - M1a_lnL)
    if LRT1 < 0:
        LRT1 = False
    
    LRT2 = 2*(M8_lnL - M7_lnL)
    if LRT2 < 0:
        LRT2 = False
    
    # run chi2 statistics to get p values
    if LRT1 != False:
        M2a_M1a = cdf_chi2(df_M2a_M1a, LRT1)
    else:
        LRT1 = None
        M2a_M1a = None
        message = "\n M1a model is more likely than M2a model."
        print(message)
        logging.info(message)
    
    if LRT2 != False:
        M8_M7 = cdf_chi2(df_M8_M7, LRT2)
    else:
        LRT2 = None
        M8_M7 = None
        message = "\n M7 model is more likely than M8 model."
        print(message)
        logging.info(message)
    
    message = "\n PAML LRTs for protein %s are : M2a v M1a - %s and M8 v M7 - %s" \
            % (protein, str(LRT1), str(LRT2))
    print(message)
    logging.info(message)
    
    # restore working directory
    os.chdir(wdir)
    
    return M2a_M1a, M8_M7
    

def run_RAxML(protein, protein_folder_path, translated_path):
    tree_name =  protein + "_Tree"
    RAxML_cline = ["raxmlHPC", "-f", "a", "-s", translated_path, "-n", tree_name, \
                   "-m", "PROTGAMMAAUTO", "-p", "1985", "-x", "2020", "-#", "100"]
    result = subprocess.call(RAxML_cline)
    all_tree_files = glob.glob("RAxML*")
    # move all the RAxML files to protein folder
    for file in all_tree_files:
        file_path = os.path.abspath(file)
        shutil.move(file_path, protein_folder_path)
    
    if result == 0:
        message = "\n Best gene tree was found for protein : %s " % protein
        print(message)
        logging.info(message)
        
    else:
        message = "\n PROBLEM with making gene tree for protein : %s " % protein
        print(message)
        logging.info(message)
    
    best_tree_path = protein_folder_path + "/RAxML_bestTree." + protein + "_Tree"
    
    # remove dashes from the tree to match species names in phylip file
    with open(best_tree_path, "r") as f:
        tree = f.read()   
        tree_no_dashes = tree.replace("_", "")
    
    with open(best_tree_path, "w") as f:
        f.write(tree_no_dashes)
    
    return best_tree_path

def get_original_cds(wdir, protein_name, original_species):

    # open according cds fasta file
    with open(wdir + "Coding_sequences/" + protein_name + "_" + original_species + "_cds.fasta", "r") as f:
        sequence = ""
        cds = f.readlines()
        for line in cds[1:]:
            sequence = sequence + line.rstrip("\n")
    return sequence

def translated_frameshift_checkpoint(wdir, seqs, protein, original_species):
    from Bio import pairwise2
    
    # take the cds of original species post-Gblocks (always the first one)
    original_post_Gblocks = seqs[0][1]
    # make a Seq object
    original_post_Gblocks_seq = Seq(original_post_Gblocks)
    # translate the Seq object
    original_post_Gblocks_seq_object = original_post_Gblocks_seq.translate()
    # make SeqRecord
    record_post_Gblocks = SeqRecord(original_post_Gblocks_seq_object)

    # take the cds coming from ensembl
    original_cds = get_original_cds(wdir, protein, original_species)
    # make a Seq object
    original_cds_seq = Seq(original_cds)
    # translate the Seq object
    original_cds_seq_object = original_cds_seq.translate()
    # make SeqRecord
    original_cds_record = SeqRecord(original_cds_seq_object)

    # perform pairwise alignment
    aln = pairwise2.align.globalxx(record_post_Gblocks.seq, \
                                        original_cds_record.seq)
    
    # make a dict storing the alignments
    d = {}
    for i in range(len(aln[0][0])):
        d[i] = ()
    
    # fill the dict with paired positions
    for position in d:
        d[position] = (aln[0][0][position], aln[0][1][position])
    
    # flag any non-synonymous substitutions that indicate frameshifts
    frameshift_positions = {}
    translated_frameshift = False
    for position, pair in d.items():
        if pair[0] != pair[1] and pair[0] != "-" and pair[1] != "-":
            # sequence index starts at 1 so need to "+1" from the python indexing
            frameshift_positions[position + 1] = pair
            translated_frameshift = True
        
    return translated_frameshift, frameshift_positions
    

def translate_Gblocks(wdir, protein_folder_path, out_Gblocks, protein, original_species):
    
    filepath_to_translate = protein_folder_path + "/" + out_Gblocks
    # read the fasta alignment file 
    alignment = AlignIO.read(filepath_to_translate, "fasta")
    # make a list of all headers and sequences as strings
    seqs = [fasta_reader.read_fasta_record(record.format("fasta")) for record in alignment]
    # make a filepath
    translated_path = filepath_to_translate.rstrip(".fasta") + "_translated.fasta"
    
    # check frameshifts in the original species aa seq post-Gblocks
    
    
    # doesnt work if pairwise2 doesnt return alignement (Clasp1)
    
    translated_frameshift, frameshift_positions = translated_frameshift_checkpoint(wdir, seqs, protein, original_species)
    if translated_frameshift == True:
        side_note = " WARNING : Frameshift in original cds post-Gblocks detected: \n%s" % frameshift_positions
        print(side_note)
        logging.info(side_note)
    
    with open(translated_path, "w") as f:
        # translate each sequence
        for s in seqs:            
            sequence = s[1]
            # make a Seq object
            coding_dna = Seq(sequence)
            # translate the Seq object
            seq_object = coding_dna.translate()
            # make SeqRecord
            record = SeqRecord(seq_object)
            # get a translated string using format method (only for SeqRecords)
            translated = record.format("fasta").lstrip("><unknown id> <unknown description\n")
            # write headers and sequences
            f.write(s[0] + "\n")
            f.write(translated)
    # return the path to the translated alignment
    return translated_path

def run_seqret(protein, protein_folder_path, out_Gblocks):
    in_filepath = protein_folder_path + "/" + out_Gblocks
    phylip_path = in_filepath + ".phy"
    seqret_cline = ["seqret", "-sequence", in_filepath, \
                    "-osformat2", "phylipnon", "-outseq", phylip_path]
    result = subprocess.call(seqret_cline)
    if result == 0:
        message = "\n Phylip format was created for protein : %s " % protein
        print(message)
        logging.info(message)
        return phylip_path
    else:
        message = "\n PROBLEM with making phylip for protein : %s " % protein
        print(message)
        logging.info(message)
      

def STOP_remover(protein_folder_path, no_dashes_out_MAFFT, protein_name):
    
    to_trim = {}
    max_trim_position = 30
    c_term_bp_to_trim = 0
    path_to_trim = protein_folder_path + "/" + no_dashes_out_MAFFT
    with open(path_to_trim, "r") as f:
        file = f.readlines()
        
        for line in file:
            
            if line.startswith(">"):
                seq = ""
                header = line.rstrip("\n")
                to_trim[header] = seq
            
            else:
                s = line.rstrip("\n")
                seq = seq + s
                to_trim[header] = seq
    
    for species, cds in to_trim.items():
        
        # split the end of each cds to look for early STOP codon
        c_term = cds[-max_trim_position:]
        c_term_codons = [c_term[i:i+3] for i in range(0, len(c_term), 3)]
            
        for codon in c_term_codons:
            if codon.upper() in ["TAA", "TGA", "TAG"]:
                # find how many bp to remove from c_term
                new_c_term_bp_to_trim = (len(c_term_codons) - c_term_codons.index(codon)) * 3
                
                # make sure that the earliest STOP is the final trim position
                if new_c_term_bp_to_trim > c_term_bp_to_trim:
                    c_term_bp_to_trim = new_c_term_bp_to_trim
    
    final_cds_no_STOP = {}
    for species, cds in to_trim.items():
        # subtract number of dashes to remove longer more bases
        final_cds_no_STOP[species] = cds[:-c_term_bp_to_trim]
        

    final_cds_file_no_STOP = "aligned_MAFFT_" + protein_name + "_final_no_STOP.fasta"
    with open(final_cds_file_no_STOP, "w") as f:
    
        for species, cds in final_cds_no_STOP.items():
            f.write(species + "\n")
            f.write(cds + "\n")
    
    return final_cds_file_no_STOP


def post_Gblocks_STOP_remover(wdir, protein_folder_path, out_Gblocks, protein): 
    # checks if artificial STOP was formed after Gblocks
    
    post_Gblocks_path = protein_folder_path + "/" + out_Gblocks
    
    all_cds = {}
    all_cds_no_STOP = {}
    artificial_STOP_codon_positions = set() # doesnt allow repetitions
    species_with_STOP = []
    
    # read species and cds into a dict to ease search
    with open(post_Gblocks_path, "r") as f:
        file = f.readlines()
        
        for line in file:
            
            if line.startswith(">"):
                seq = ""
                header = line.rstrip("\n")
                all_cds[header] = seq
            
            else:
                s = (line.rstrip("\n").replace(" ", "")).upper()
                seq = seq + s
                all_cds[header] = seq
    
    # chop cds into codons (all are in frame) and mark STOPs
    for header, cds in all_cds.items():

        codons_list = [cds[i:i+3] for i in range(0, len(cds), 3)]
        STOP_positions = [position for position, codon in enumerate(codons_list) if codon in ["TAA", "TGA", "TAG"]]

        if STOP_positions:
            
            species = header.replace(">", "")
            species_with_STOP.append(species)
            
            for position in STOP_positions:
                artificial_STOP_codon_positions.add(position)
    
    # edit all coding sequences removing STOP columns
    for header, cds in all_cds.items():
        
        all_codons = [cds[i:i+3] for i in range(0, len(cds), 3)]
        cds_no_STOP = ""
        
        # join all the codon except STOP positions
        all_cds_no_STOP[header] = "".join([cds_no_STOP + codon for position, codon in enumerate(all_codons) \
                                   if position not in artificial_STOP_codon_positions])
    
    # overwrite the previous out_Gblocks file
    with open(post_Gblocks_path, "w") as f:
        
        for header, cds in all_cds_no_STOP.items():
            f.write(header + "\n")
            f.write(cds + "\n")
    
    # issue a warning if STOPs were present
    if species_with_STOP:
        
        STOP_positions = sorted([(position * 3) + 1 for position in artificial_STOP_codon_positions])
        
        message = "\n !!! WARNING !!! - CDS for: %s contains early STOP codons "\
        "starting at nucleotide position in Gblocks_final_no_STOP alignment: %s in species: %s \n           "\
        "----> they were removed in each species forcing conserved alignment" \
                % (protein, STOP_positions, species_with_STOP) 
        print(message)
        logging.info(message)
    
    return post_Gblocks_path


def run_Gblocks(final_cds_file_no_STOP, protein, result_path):
    in_filepath = result_path + protein + "/" + final_cds_file_no_STOP
    out_Gblocks = "aligned_MAFFT_Gblocks_" + protein + "_final_no_STOP.fasta"
    # run Gblocks with options: codon ("-t=c") and dont save html file ("-p=n")
    Gblocks_cline = ["Gblocks", in_filepath, "-t=c", "-p=n"]
    result = subprocess.run(Gblocks_cline, capture_output=True)
    # get stdout -> decode to string from bit -> split by lines and get the 5th and 6th
    message = (result.stdout.decode('utf-8').split("\n"))[5:7]
    # log stdout
    print(message)
    logging.info(message)
    # add fasta extension
    os.rename(in_filepath + "-gb", out_Gblocks)
    
    # returns the filename after Gblocks
    return out_Gblocks

              
def align_final_cds(protein, final_cds_file, result_path):
    
    in_filepath = result_path + protein + "/" + protein + "_final.fasta"
    out_MAFFT = "aligned_MAFFT_" + final_cds_file
    # run mafft
    mafft_cline = MafftCommandline(input=in_filepath)
    # record standard output and standard error
    stdout, stderr = mafft_cline()
    # make a post-MSA file using out_filename
    with open(out_MAFFT, "w") as f:
        f.write(stdout)
    
    # returns the filename after MAFFT
    return out_MAFFT


"""

def remove_non_ACGT_codons(result_path, out_Gblocks):
    
    non_ACGT = ["N","Y","R","W","S","K","M","D","H","V","B","X"]
    cleaned_Gblocks = open(out_Gblocks[:-6] + "_cleaned.fasta") 
    out_Gblocks = result_path + "aligned_MAFFT_Gblocks_Terf1_final.fasta"
    with open(out_Gblocks, "r") as f:
        line_nr = 0
        file = f.readlines()
        for line in file:
            line_nr += 1
            if line.startswith('>'):
                pass
                #cleaned_Gblocks.write(line)
            elif any(n in line.upper() for n in non_ACGT):
                print(line)
                print(line_nr)
            else:
                pass
                #cleaned_Gblocks.write(line)
"""



