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

Runs blast using NCBI tblastn.

"""

from freeda import genomes_preprocessing, pyinstaller_compatibility
import subprocess
import os
import glob
import time
import logging


def run_blast(wdir, ref_species, all_genes, final_excluded_species=None):
    """Runs tblastn based on NCBI makedatabase routine."""

    database_path = wdir + "Genomes/"
    query_path = wdir + "Blast_input/"
    output_path = wdir + "Blast_output/"
    form = "6 qseqid means sseqid means qstart means qend means sstart means send means evalue means " \
           "bitscore length means pident means mismatch means gapopen means qlen means slen means"

    all_genomes = genomes_preprocessing.get_names(wdir, ref_species, final_excluded_species)
    genomes = [names[1] for names in all_genomes]

    # clear Blast_output folder
    all_old_blast_output_files = glob.glob(os.path.join(output_path, "*.txt"))
    for file in all_old_blast_output_files:
        os.remove(file)

    # make sure database for each genome already exists or is successfully built
    for genome in genomes:
        genome_file_database = check_genome_present(wdir, ref_species, database_path, genome, ref_genome=False)
        # failed to build database
        if genome_file_database is False:
            return None

    # PYINSTALLER: Add path to os path variable.
    tblastn_path = pyinstaller_compatibility.resource_path("tblastn")

    # perform blast
    for genome in genomes:
        for gene in all_genes:
            database = database_path + genome + ".fasta"
            query = query_path + gene + "_" + ref_species + "_protein.fasta"
            output = output_path + gene + "_" + genome + ".txt"
            to_blast = [tblastn_path, "-db", database, "-query", query,
                        "-out", output, "-outfmt", form, "-num_threads", "8"]
            message = "\nPerforming tblastn for gene: %s from genome: %s\n" % (gene, genome)
            logging.info(message)
            subprocess.call(to_blast)

    print("\ntblastn txt files have been generated.")
    
    return


def remove_empty_files(database_path):
    """Checks if there are any empty files (failed downloads) from previous runs and eliinates them."""

    # get all files in the genomes directory
    all_files = [database_path + f for f in os.listdir(database_path) if os.path.isfile(os.path.join(database_path, f))]
    # check is any of them is empty
    for f in all_files:
        if not os.stat(f).st_size:
            os.remove(f)


def check_genome_downloads(ref_species, database_path, genome, zip=False):
    """Removes failed downloads -> .zip and .fasta files."""

    # check for zip files (should be absent)
    if zip:
        if os.path.exists(database_path + genome + ".zip"):
            os.remove(database_path + genome + ".zip")
            return False
        else:
            return True

    # fasta files size (-100 000 bytes)
    rodents = {"MusMusculus_genome": 2762231583,  # zip 870318080
             "MusSpicilegus_genome": 2530495171,  # zip 786432000
             "MusSpretus_genome": 2658802980,  # zip 781189120
             "MusCaroli_genome": 2585244504,  # zip 754974720,
             "MusMinutoides_genome": 2876597405,  # zip 807403520
             "MusPahari_genome": 2506109361,  # zip 744488960
             "ApodemusSylvaticus_genome": 3860843764,  # zip 744488960
             "ApodemusSpeciosus_genome": 3576911245,  # zip 723517440
             "HylomyscusAlleni_genome": 2268841987,  # zip 734003200
             "PraomysDelectorum_genome": 2461069712,  # zip 786432000
             "MastomysNatalensis_genome": 2548983787,  # zip 807403520
             "MastomysCoucha_genome": 2538410655,  # zip 765460480
             "GrammomysDolichurus_genome": 2168341549,  # zip 692060160
             "GrammomysSurdaster_genome": 2444834577,  # zip 786432000
             "ArvicanthisNiloticus_genome": 2528074480,  # zip 807403520
             "RhabdomysDilectus_genome": 2288892761,  # zip 713031680
             "RhynchomysSoricoides_genome" : 2211968598,  # zip 713031680
             "RattusRattus_genome": 2412685692,  # zip 786432000
             "RattusNorvegicus_genome": 2906067117}  # zip 891289600

    primates = {"HomoSapiens_genome": 3138414431,
             "PanTroglodytes_genome": 3062328323,
             "GorillaGorilla_genome": 3083663602,
             "PongoAbelli_genome": 3103792067,
             "NomascusLeucogenys_genome": 2879637870,
             "HylobatesMoloch_genome": 2885571806,
             "SymphalangusSyndactylus_genome": 2808876019,  # added 08/17/22
             "HoolockLeuconedys_genome": 2820095160,  # added 08/17/22
             "CercopithecusMona_genome": 2939174695,
             "MacacaMulatta_genome": 3075206655,
             "PapioAnubis_genome": 2906552752,
             "ChlorocebusSabaeus_genome": 2975174619,
             "ErythrocebusPatas_genome": 3136127270,  # added 08/17/22
             "MandrillusSphinx_genome": 2872009988,  # added 08/17/22
             "LophocebusAterrimus_genome": 2940078738,  # added 08/17/22
             "TrachypithecusFrancoisi_genome": 2938846753,
             "PiliocolobusTephrosceles_genome": 3080829369,
             "RhinopithecusStrykeri_genome": 2971349299,  # added 08/17/22
             "PygathrixNigripes_genome": 2932771300,  # added 08/17/22
             "PitheciaPithecia_genome": 3060979898,
             "AotusNancymaae_genome": 2900243235,
             "PlecturocebusDonacophilus_genome": 3112365802,
             "AlouattaPalliata_genome": 3188439879,
             "CallithrixJacchus_genome": 2761618951,
             "SaimiriBoliviensis_genome": 2684668098,
             "AtelesGeoffroyi_genome": 3020749703,
             "CebusAlbifrons_genome": 2910572252,  # added 08/17/22
             "SapajusApella_genome": 2794074315}  # added 08/17/22

    carnivora = {"CanisFamiliaris_genome": 2705909865,
           "SpeothosVenaticus_genome": 2351270082,
           "VulpesFerrilata_genome": 2409426090,
           "UrsusMaritimus_genome": 2359879905,
           "UrsusAmericanus_genome": 2631398004,
           "MiroungaLeonina_genome": 2447526042,
           "OdobenusRosmarus_genome": 2430488797,
           "ZalophusCalifornianus_genome": 2439711162,
           "AilurusFulgens_genome": 2373290453,
           "ProcyonLotor_genome": 2283274334,
           "MelesMeles_genome": 2772983809,
           "GuloGulo_genome": 2287813027,
           "MustelaNigripes_genome": 2418506728,
           "LutraLutra_genome": 2365629721,
           "SpilogaleGracilis_genome": 2540140529,
           "ParadoxurusHermaphroditus_genome": 2559496470,
           "CryptoproctaFerox_genome": 2539739378,
           "SuricataSuricatta_genome": 2389361575,
           "HyaenaHyaena_genome": 2512114800,
           "PantheraTigris_genome": 2437681936,
           "LynxRufus_genome": 2469754470,
           "FelisCatus_genome": 2486898359}

    phasianidae = {"GallusGallus_genome": 1248683961,
               "BambusicolaThoracicus_genome": 1061115573,
               "PavoCristatus_genome": 0,  # added 08/17/22
               "PavoMuticus_genome": 1074499003,
               "MeleagrisGallopavo_genome": 1158971686,
               "CentrocercusUrophasianus_genome": 1023355947,
               "CentrocercusMinimus_genome": 0,  # added 08/17/22
               "LyrurusTetrix_genome": 748849483,
               "LagopusLeucura_genome": 1029379726,
               "LagopusMuta_genome": 0,  # added 08/17/22
               "TympanuchusCupido_genome": 997281378,
               "ChrysolophusPictus_genome": 1038257264,
               "PhasianusColchicus_genome": 1034344338,
               "CrossoptilonMantchuricum_genome": 1026225886,
               "SyrmaticusMikado_genome": 1086380439,
               "LophuraNycthemera_genome": 0,  # added 08/17/22
               "AlectorisRufa_genome": 1041225438,
               "CoturnixJaponica_genome": 939336855}

    # define clade
    if ref_species in ["Mm", "Rn"]:
        genomes = rodents
    elif ref_species in ["Hs"]:
        genomes = primates
    elif ref_species in ["Cf", "Fc"]:
        genomes = carnivora
    elif ref_species in ["Gg"]:
        genomes = phasianidae

    genome_file = False
    # check if genome fasta file was downloaded and unpacked successfully
    for root, dirs, files in os.walk(database_path, topdown=False):
        for file in files:
            if file == genome + ".fasta":
                # check size
                if os.stat(database_path + genome + ".fasta").st_size < genomes[genome]:
                    # partial fasta file -> remove
                    message = "\n...NOTE... : Partial genome file : %s.fasta detected: " \
                              "\n     -> removing..." % genome
                    logging.info(message)
                    os.remove(database_path + genome + ".fasta")
                    return genome_file
                # complete fasta file
                else:
                    genome_file = True
                    return genome_file

    # no fasta file
    return genome_file


def check_genome_present(wdir, ref_species, database_path, genome, ref_genome=False):
    """Checks if a given genome is present. Unpacks and unzips genomes downloaded from NCBI Assembly.
    Non-ncbi assemblies must be prepared as ".fasta" files conform with "genomes.txt" names.
    It also looks for reference genome if key-only argument reference_genome is invoked."""

    # check for empty files (failed downloads)
    remove_empty_files(database_path)

    # if function called by input extractor module not tblastn module
    if ref_genome is True:
        ref_genome_path = database_path

    # get info on all files available
    all_files = []
    for root, dirs, files in os.walk(database_path, topdown=False):
        for f in files:
            all_files.append(f)

    ### --- Getting the non-reference genomes --- ###

    if not ref_genome:

        genome_file_databases = False

        # check if genome database is present (most cases)
        if genome + ".fasta.nin" in all_files:
            message = "\nGenome : %s blast database already exists" % genome
            logging.info(message)
            genome_file_databases = True
            return genome_file_databases

        # check if a compressed zip file from previous run is present (and remove)
        if not check_genome_downloads(ref_species, database_path, genome, zip=True):
            message = "\n...NOTE... : Detected compressed genome file : %s.zip from a previous run" \
                      "\n   -> removed..." % genome
            logging.info(message)

        # complete fasta file present but no databases
        if check_genome_downloads(ref_species, database_path, genome):
            message = "\n...NOTE... : Genome : %s.fasta file already exists but blast databases are missing" % genome
            logging.info(message)

            # proceed to making databases
            make_blast_database(database_path, genome)

            # validate that the database was generated
            all_files = []
            for root, dirs, files in os.walk(database_path, topdown=False):
                for f in files:
                    all_files.append(f)

            # failed to build databases
            if genome + ".fasta.nin" not in all_files:
                message = "\n...FATAL_ERROR... : Genome : %s database failed to build" \
                              "\n     -> exiting the pipeline now..." % genome
                logging.info(message)
                return genome_file_databases

            # successful build
            else:
                message = "\nGenome : %s was downloaded and decompressed successfully" \
                              "\n        -> blast databases were built successfully" % genome
                logging.info(message)
                genome_file_databases = True
                return genome_file_databases

        # databases, zip and fasta files are all absent
        else:
            message = "\n...NOTE... : Did not detect genome or blast databases for : %s" \
              "\n -> downloading and decompressing the genome (it might take a couple of minutes)...\n" % genome
            logging.info(message)

            # download genome
            all_genomes = genomes_preprocessing.get_names(wdir, ref_species)
            accession_nr = [names[2] for names in all_genomes if genome in names][0]
            download_genome(genome, accession_nr, database_path)

            # check if zip file is still present (should not be)
            if not check_genome_downloads(ref_species, database_path, genome, zip=True):
                message = "\n...FATAL_ERROR... : Genome : %s failed to download or decompress" \
                          "   \n        -> likely Internet connection was disrupted" \
                          "\n               -> please run FREEDA again..." % genome
                logging.info(message)
                return genome_file_databases

            # check if fasta file is still partial/empty (should not be)
            if not check_genome_downloads(ref_species, database_path, genome):
                message = "\n...FATAL_ERROR... : Partial or empty genome file : %s.fasta" \
                          "   \n        -> likely Internet connection was disrupted" \
                          "\n               -> please run FREEDA again..." % genome
                logging.info(message)
                return genome_file_databases

            # proceed to making databases
            make_blast_database(database_path, genome)

            # validate that the database was generated
            all_files = []
            for root, dirs, files in os.walk(database_path, topdown=False):
                for f in files:
                    all_files.append(f)

            # failed to build databases
            if genome + ".fasta.nin" not in all_files:
                message = "\n...FATAL_ERROR... : Genome : %s database failed to build" \
                      "\n     -> exiting the pipeline now..." % genome
                logging.info(message)
                return genome_file_databases

            # successful build
            else:
                message = "\nGenome : %s was downloaded and decompressed successfully" \
                      "\n        -> blast databases were built successfully" % genome
                logging.info(message)
                genome_file_databases = True
                return genome_file_databases

    ### ---- Getting the reference genome ---- ###

    if ref_genome:

        ref_genome_file = False

        # check if ref genome present (most of the cases)
        if check_genome_downloads(ref_species, database_path, genome):
            message = "\nReference genome : %s is present" % genome
            logging.info(message)
            ref_genome_file = True
            return ref_genome_file

        # check if a compressed zip file from previous run is present (and remove)
        if not check_genome_downloads(ref_species, database_path, genome, zip=True):
            message = "\n...NOTE... : Detected compressed reference genome file : %s.zip from a previous run" \
                        "\n      -> removed..." % genome
            logging.info(message)

        # no complete fasta file detected
        if not check_genome_downloads(ref_species, database_path, genome):
            message = "\n...NOTE... : Did not detect a complete reference genome file : %s.fasta" \
                    "\n    -> downloading and decompressing the genome (it might take a couple of minutes)...\n" % genome
            logging.info(message)

            # download genome
            all_genomes = genomes_preprocessing.get_names(wdir, ref_species, ref_genome=True)
            accession_nr = all_genomes[2]
            download_genome(genome, accession_nr, ref_genome_path)

            # check if fasta size is correct:
            if check_genome_downloads(ref_species, database_path, genome):
                message = "\nReference genome : %s was downloaded and decompressed successfully" % genome
                logging.info(message)
                ref_genome_file = True
                return ref_genome_file

            # failed to download or decompress the ref genome
            else:
                message = "\n...FATAL_ERROR... : Reference genome : %s downloading or decompression failed" \
                            "   \n      -> likely Internet connection was disrupted" \
                            "\n             -> please run FREEDA again..." % genome
                logging.info(message)
                return ref_genome_file

    else:
        message = "\nElse statement when making blast database\n"
        logging.info(message)
        print("Else statement when making blast database")


def download_genome(genome, accession_nr, database_path):
    """Downloads and unzips genomes using NCBI Datasets CLI, makes one fasta file"""

    # download genome
    start_time = time.time()
    filepath_1 = database_path + genome + ".zip"
    # need to exclude genomic cds because its also ".fna" which confuses the concatenation
    cmd1 = [pyinstaller_compatibility.resource_path("datasets"),
            "download", "genome", "accession", accession_nr, "--exclude-genomic-cds", "--filename", filepath_1]
    subprocess.call(cmd1, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))  # mute output and cautions

    # unzip all chromosomes into single file
    filepath_2 = database_path + genome + ".fasta"
    cmd2 = [pyinstaller_compatibility.resource_path("unzip"), "-pq", filepath_1, "cat", "*/*.fna", filepath_2]
    with open(filepath_2, "w") as outfile:
        subprocess.call(cmd2, stdout=outfile, stderr=open(os.devnull, 'wb'))  # redirect to file and mute cautions

    stop_time = time.time()
    message = "         -> Done : in %s min" % ((stop_time - start_time) / 60)
    logging.info(message)

    os.remove(filepath_1)


def make_blast_database(database_path, genome):
    """ Makes a blast database based on genome fasta file."""

    genome_file = genome + ".fasta"
    make_database_nucl = [pyinstaller_compatibility.resource_path("makeblastdb"),
                          "-in", database_path + genome_file, "-dbtype", "nucl"]
    make_database_prot = [pyinstaller_compatibility.resource_path("makeblastdb"),
                          "-in", database_path + genome_file, "-dbtype", "prot"]

    message = "\n                  Building blast database for genome : %s ..." % genome
    logging.info(message)
    subprocess.call(make_database_nucl, stdout=open(os.devnull, 'wb'))  # mute log info
    subprocess.call(make_database_prot, stdout=open(os.devnull, 'wb'))  # mute log info
