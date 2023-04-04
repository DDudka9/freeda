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
import platform


def run_blast(wdir, ref_species, all_genes, final_excluded_species=None, subgroup=None):
    """Runs tblastn based on NCBI makedatabase routine."""

    database_path = wdir + "Genomes/"
    query_path = wdir + "Blast_input/"
    output_path = wdir + "Blast_output/"
    form = "6 qseqid means sseqid means qstart means qend means sstart means send means evalue means " \
           "bitscore length means pident means mismatch means gapopen means qlen means slen means"

    # allow increased threads for iOS systems
    if platform.uname().system != "Darwin":
        threads = 1
    else:
        threads = 4

    all_genomes = genomes_preprocessing.get_names(wdir, ref_species, final_excluded_species, subgroup)
    genomes = [names[1] for names in all_genomes]

    # clear Blast_output folder
    all_old_blast_output_files = glob.glob(os.path.join(output_path, "*.txt"))
    for file in all_old_blast_output_files:
        os.remove(file)

    # make sure database for each genome already exists or is successfully built
    for genome in genomes:
        genome_file_database = check_genome_present(wdir, ref_species, database_path, genome,
                                                    final_excluded_species, subgroup, ref_genome=False)
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
                        "-out", output, "-outfmt", form, "-num_threads", str(threads)]  # "-num_threads", "4" removed (11/04/22)
            message = "\nPerforming tblastn for gene: %s from genome: %s\n" % (gene, genome)
            logging.info(message)
            subprocess.call(to_blast)

    print("\ntblastn txt files have been generated.")
    
    return


def remove_empty_files(database_path):
    """Checks if there are any empty files (failed downloads) from previous runs and eliminates them."""

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
             'PanTroglodytes_genome.fasta': 3062328323,
             'PanPaniscus_genome.fasta': 3312469357,
             'GorillaGorilla_genome.fasta': 3083663602,
             'PongoAbelli_genome.fasta': 3103792067,
             'PongoPygmaeus_genome.fasta': 3290024527,
             'NomascusLeucogenys_genome.fasta': 2879637870,
             'NomascusSiki_genome.fasta': 2818036002,
             'HylobatesMoloch_genome.fasta': 2885571806,
             'HylobatesPileatus_genome.fasta': 2887024203,
             'SymphalangusSyndactylus_genome.fasta': 3222691930,
             'HoolockLeuconedys_genome.fasta': 2819995160,
             'CercopithecusMona_genome.fasta': 2939174695,
             'MacacaMulatta_genome.fasta': 3075206655,
             'PapioAnubis_genome.fasta': 2906552752,
             'ChlorocebusSabaeus_genome.fasta': 2975174619,
             'TrachypithecusFrancoisi_genome.fasta': 2938846753,
             'PiliocolobusTephrosceles_genome.fasta': 3080829369,
             'PitheciaPithecia_genome.fasta': 3060979898,
             'AotusNancymaae_genome.fasta': 2900243235,
             'PlecturocebusDonacophilus_genome.fasta': 3112365802,
             'AlouattaPalliata_genome.fasta': 3188439879,
             'CallithrixJacchus_genome.fasta': 2761618951,
             'SaimiriBoliviensis_genome.fasta': 2684668098,
             'AtelesGeoffroyi_genome.fasta': 3020749703}

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
           "MungosMungo_genome": 0,
           "HyaenaHyaena_genome": 2512114800,
           "PantheraTigris_genome": 2437681936,
           "NeofelisNebulosa_genome": 0,
           "PumaConcolor_genome": 0,
           "AcinonyxJubatus_genome": 0,
           "LynxRufus_genome": 2469754470,
           "PrionailurusViverrinus_genome": 0,
           "OtocolobusManul_genome": 0,
           "FelisCatus_genome": 2486898359}

    phasianidae = {"GallusGallus_genome": 1248683961,
               "BambusicolaThoracicus_genome": 1061115573,
               "PavoCristatus_genome": 1059867822,
               "PavoMuticus_genome": 1074499003,
               "MeleagrisGallopavo_genome": 1158971686,
               "CentrocercusUrophasianus_genome": 1023355947,
               "CentrocercusMinimus_genome": 1012645442,
               "LyrurusTetrix_genome": 748849483,
               "LagopusLeucura_genome": 1029379726,
               "LagopusMuta_genome": 1039604974,
               "TympanuchusCupido_genome": 997281378,
               "ChrysolophusPictus_genome": 1038257264,
               "PhasianusColchicus_genome": 1034344338,
               "CrossoptilonMantchuricum_genome": 1026225886,
               "SyrmaticusMikado_genome": 1086380439,
               "LophuraNycthemera_genome": 1027348961,
               "AlectorisRufa_genome": 1041225438,
               "CoturnixJaponica_genome": 939336855}

    flies = {"DrosophilaMelanogaster_genome": 145550290,
             "DrosophilaSimulans_genome": 135164843,
             "DrosophilaMauritiana_genome": 131710764,
             "DrosophilaSechellia_genome": 149708744,
             "DrosophilaYakuba_genome": 153141101,
             "DrosophilaSantomea_genome": 148564044,
             "DrosophilaTeissieri_genome": 146510489,
             "DrosophilaOrena_genome": 185126277,
             "DrosophilaErecta_genome": 137303873,
             "DrosophilaEugracilis_genome": 167000171,
             "DrosophilaSubpulchrella_genome": 271471590,
             "DrosophilaSuzukii_genome": 271313394,
             "DrosophilaBiarmipes_genome": 187564475,
             "DrosophilaTakahashii_genome": 167510332,
             "DrosophilaFicusphila_genome": 169919136,
             "DrosophilaCarrolli_genome": 234041655,
             "DrosophilaRhopaloa_genome": 195850575,
             "DrosophilaKurseongensis_genome": 208948327,
             "DrosophilaFuyamai_genome": 231955431,
             "DrosophilaElegans_genome": 180649414}

    # define clade
    if ref_species in ["Mm", "Rn"]:
        genomes = rodents
    elif ref_species in ["Hs"]:
        genomes = primates
    elif ref_species in ["Cf", "Fc"]:
        genomes = carnivora
    elif ref_species in ["Gg"]:
        genomes = phasianidae
    elif ref_species in ["Dme"]:
        genomes = flies

    genome_file = False

    # removed the e

    # check if genome fasta file was downloaded and unpacked successfully
    for root, dirs, files in os.walk(database_path, topdown=False):
        for file in files:
            if file == genome + ".fasta":

                # THIS WAS ADDED TEMPORARILY 02/16/2023
                #message = "GENOME file: %s present of size: %s " % (file,
                #                            str(os.stat(database_path + genome + ".fasta").st_size))
                #logging.info(message)
                #genome_file = True
                #return genome_file

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


def check_genome_present(wdir, ref_species, database_path, genome, final_excluded_species=None,
                         subgroup=None, ref_genome=False):
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
        if check_blast_database(genome, database_path):
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
            #all_files = []
            #for root, dirs, files in os.walk(database_path, topdown=False):
            #    for f in files:
            #        all_files.append(f)

            # failed to build databases
            if not check_blast_database(genome, database_path):
                genome_file_databases = False
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
            all_genomes = genomes_preprocessing.get_names(wdir, ref_species, final_excluded_species, subgroup)
            accession_nr = [names[2] for names in all_genomes if genome in names][0]
            download_genome(genome, accession_nr, database_path)

            # check if zip file is still present (should not be)
            if not check_genome_downloads(ref_species, database_path, genome, zip=True):
                message = "\n...FATAL ERROR... : Genome : %s failed to download or decompress" \
                          "   \n        -> likely Internet connection was disrupted" \
                          "\n               -> please run FREEDA again..." % genome
                logging.info(message)
                return genome_file_databases

            # check if fasta file is still partial/empty (should not be)
            if not check_genome_downloads(ref_species, database_path, genome):
                message = "\n...FATAL ERROR... : Partial or empty genome file : %s.fasta" \
                          "   \n        -> likely Internet connection was disrupted" \
                          "\n               -> please run FREEDA again..." % genome
                logging.info(message)
                return genome_file_databases

            # proceed to making databases
            make_blast_database(database_path, genome)

            # validate that the database was generated
            #all_files = []
            #for root, dirs, files in os.walk(database_path, topdown=False):
            #    for f in files:
            #        all_files.append(f)

            # failed to build databases
            if not check_blast_database(genome, database_path):
                genome_file_databases = False
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
                message = "\n...FATAL ERROR... : Reference genome : %s downloading or decompression failed" \
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
    #make_database_prot = [pyinstaller_compatibility.resource_path("makeblastdb"),  # removed 11/04/22
    #                      "-in", database_path + genome_file, "-dbtype", "prot"]

    message = "\n                  Building blast database for genome : %s ..." % genome
    logging.info(message)
    subprocess.call(make_database_nucl, stdout=open(os.devnull, 'wb'))  # mute log info
    #subprocess.call(make_database_prot, stdout=open(os.devnull, 'wb'))  # mute log info  # removed 11/04/22


def check_blast_database(genome, database_path):
    """Checks size of all expected blast databae files. Files and sizes based on BLAST 2.12.0 version."""

    #  -100 bytes from each expected size
    possible_databases = {
        'AilurusFulgens': [('AilurusFulgens_genome.fasta.ntf', 16284),
                            ('AilurusFulgens_genome.fasta.nto', 46260),
                            ('AilurusFulgens_genome.fasta.not', 138976),
                            ('AilurusFulgens_genome.fasta.ndb', 20380),
                            ('AilurusFulgens_genome.fasta.nsq', 586879275),
                            ('AilurusFulgens_genome.fasta.nhr', 1784478),
                            ('AilurusFulgens_genome.fasta.nin', 139156)],
         'AlectorisRufa': [('AlectorisRufa_genome.fasta.ndb', 20380),
                           ('AlectorisRufa_genome.fasta.ntf', 16284),
                           ('AlectorisRufa_genome.fasta.nhr', 1653017),
                           ('AlectorisRufa_genome.fasta.nsq', 257072701),
                           ('AlectorisRufa_genome.fasta.not', 127084),
                           ('AlectorisRufa_genome.fasta.nin', 127256),
                           ('AlectorisRufa_genome.fasta.nto', 42296)],
         'AlouattaPalliata': [('AlouattaPalliata_genome.fasta.ndb', 20380),
                              ('AlouattaPalliata_genome.fasta.ntf', 16284),
                              ('AlouattaPalliata_genome.fasta.nhr', 188084951),
                              ('AlouattaPalliata_genome.fasta.nsq', 759308424),
                              ('AlouattaPalliata_genome.fasta.not', 13638688),
                              ('AlouattaPalliata_genome.fasta.nin', 13638868),
                              ('AlouattaPalliata_genome.fasta.nto', 4546164)],
         'AotusNancymaae': [('AotusNancymaae_genome.fasta.ndb', 20380),
                            ('AotusNancymaae_genome.fasta.ntf', 16284),
                            ('AotusNancymaae_genome.fasta.nhr', 4682711),
                            ('AotusNancymaae_genome.fasta.nsq', 716294716),
                            ('AotusNancymaae_genome.fasta.not', 346960),
                            ('AotusNancymaae_genome.fasta.nin', 347140),
                            ('AotusNancymaae_genome.fasta.nto', 115588)],
         'ApodemusSpeciosus': [('ApodemusSpeciosus_genome.fasta.ndb', 20380),
                               ('ApodemusSpeciosus_genome.fasta.ntf', 16284),
                               ('ApodemusSpeciosus_genome.fasta.nhr', 60022095),
                               ('ApodemusSpeciosus_genome.fasta.nsq', 906584731),
                               ('ApodemusSpeciosus_genome.fasta.not', 4033396),
                               ('ApodemusSpeciosus_genome.fasta.nin', 4033576),
                               ('ApodemusSpeciosus_genome.fasta.nto', 1344400)],
         'ApodemusSylvaticus': [('ApodemusSylvaticus_genome.fasta.ndb', 20380),
                                ('ApodemusSylvaticus_genome.fasta.ntf', 16284),
                                ('ApodemusSylvaticus_genome.fasta.nhr', 90779916),
                                ('ApodemusSylvaticus_genome.fasta.nsq', 965931852),
                                ('ApodemusSylvaticus_genome.fasta.not', 6715456),
                                ('ApodemusSylvaticus_genome.fasta.nin', 6715636),
                                ('ApodemusSylvaticus_genome.fasta.nto', 2238420)],
         'ArvicanthisNiloticus': [('ArvicanthisNiloticus_genome.fasta.ndb', 20380),
                                  ('ArvicanthisNiloticus_genome.fasta.ntf', 16284),
                                  ('ArvicanthisNiloticus_genome.fasta.nhr', 308549),
                                  ('ArvicanthisNiloticus_genome.fasta.nsq', 624247348),
                                  ('ArvicanthisNiloticus_genome.fasta.not', 19060),
                                  ('ArvicanthisNiloticus_genome.fasta.nin', 19248),
                                  ('ArvicanthisNiloticus_genome.fasta.nto', 6288)],
         'AtelesGeoffroyi': [('AtelesGeoffroyi_genome.fasta.ndb', 20380),
                             ('AtelesGeoffroyi_genome.fasta.ntf', 16284),
                             ('AtelesGeoffroyi_genome.fasta.nhr', 141941982),
                             ('AtelesGeoffroyi_genome.fasta.nsq', 724918041),
                             ('AtelesGeoffroyi_genome.fasta.not', 10428808),
                             ('AtelesGeoffroyi_genome.fasta.nin', 10428988),
                             ('AtelesGeoffroyi_genome.fasta.nto', 3476204)],
         'BambusicolaThoracicus': [('BambusicolaThoracicus_genome.fasta.ndb', 20380),
                                   ('BambusicolaThoracicus_genome.fasta.ntf', 16284),
                                   ('BambusicolaThoracicus_genome.fasta.nhr', 26219322),
                                   ('BambusicolaThoracicus_genome.fasta.nsq', 258579913),
                                   ('BambusicolaThoracicus_genome.fasta.not', 1964884),
                                   ('BambusicolaThoracicus_genome.fasta.nin', 1965072),
                                   ('BambusicolaThoracicus_genome.fasta.nto', 654896)],
         'CallithrixJacchus': [('CallithrixJacchus_genome.fasta.ndb', 20380),
                               ('CallithrixJacchus_genome.fasta.ntf', 16284),
                               ('CallithrixJacchus_genome.fasta.nhr', 231337),
                               ('CallithrixJacchus_genome.fasta.nsq', 681897026),
                               ('CallithrixJacchus_genome.fasta.not', 14512),
                               ('CallithrixJacchus_genome.fasta.nin', 14692),
                               ('CallithrixJacchus_genome.fasta.nto', 4772)],
         'CentrocercusMinimus': [('CentrocercusMinimus_genome.fasta.ndb', 20380),
                                 ('CentrocercusMinimus_genome.fasta.ntf', 16284),
                                 ('CentrocercusMinimus_genome.fasta.nhr', 161933),
                                 ('CentrocercusMinimus_genome.fasta.nsq', 250011989),
                                 ('CentrocercusMinimus_genome.fasta.not', 11728),
                                 ('CentrocercusMinimus_genome.fasta.nin', 11916),
                                 ('CentrocercusMinimus_genome.fasta.nto', 3844)],
         'CentrocercusUrophasianus': [('CentrocercusUrophasianus_genome.fasta.ndb',
                                       20380),
                                      ('CentrocercusUrophasianus_genome.fasta.ntf', 16284),
                                      ('CentrocercusUrophasianus_genome.fasta.nhr', 257009),
                                      ('CentrocercusUrophasianus_genome.fasta.nsq', 252812645),
                                      ('CentrocercusUrophasianus_genome.fasta.not', 17932),
                                      ('CentrocercusUrophasianus_genome.fasta.nin', 18128),
                                      ('CentrocercusUrophasianus_genome.fasta.nto', 5912)],
         'CercopithecusMona': [('CercopithecusMona_genome.fasta.ndb', 20380),
                               ('CercopithecusMona_genome.fasta.ntf', 16284),
                               ('CercopithecusMona_genome.fasta.nhr', 300905),
                               ('CercopithecusMona_genome.fasta.nsq', 725702293),
                               ('CercopithecusMona_genome.fasta.not', 22576),
                               ('CercopithecusMona_genome.fasta.nin', 22756),
                               ('CercopithecusMona_genome.fasta.nto', 7460)],
         'ChlorocebusSabaeus': [('ChlorocebusSabaeus_genome.fasta.ndb', 20380),
                                ('ChlorocebusSabaeus_genome.fasta.ntf', 16284),
                                ('ChlorocebusSabaeus_genome.fasta.nhr', 1146289),
                                ('ChlorocebusSabaeus_genome.fasta.nsq', 734544574),
                                ('ChlorocebusSabaeus_genome.fasta.not', 82372),
                                ('ChlorocebusSabaeus_genome.fasta.nin', 82560),
                                ('ChlorocebusSabaeus_genome.fasta.nto', 27392)],
         'ChrysolophusPictus': [('ChrysolophusPictus_genome.fasta.ndb', 20380),
                                ('ChrysolophusPictus_genome.fasta.ntf', 16284),
                                ('ChrysolophusPictus_genome.fasta.nhr', 2445787),
                                ('ChrysolophusPictus_genome.fasta.nsq', 256435030),
                                ('ChrysolophusPictus_genome.fasta.not', 171436),
                                ('ChrysolophusPictus_genome.fasta.nin', 171616),
                                ('ChrysolophusPictus_genome.fasta.nto', 57080)],
         'CoturnixJaponica': [('CoturnixJaponica_genome.fasta.ndb', 20380),
                              ('CoturnixJaponica_genome.fasta.ntf', 16284),
                              ('CoturnixJaponica_genome.fasta.nhr', 324591),
                              ('CoturnixJaponica_genome.fasta.nsq', 232015025),
                              ('CoturnixJaponica_genome.fasta.not', 24040),
                              ('CoturnixJaponica_genome.fasta.nin', 24220),
                              ('CoturnixJaponica_genome.fasta.nto', 7948)],
         'CrossoptilonMantchuricum': [('CrossoptilonMantchuricum_genome.fasta.ndb',
                                       20380),
                                      ('CrossoptilonMantchuricum_genome.fasta.ntf', 16284),
                                      ('CrossoptilonMantchuricum_genome.fasta.nhr', 364678),
                                      ('CrossoptilonMantchuricum_genome.fasta.nsq', 254660092),
                                      ('CrossoptilonMantchuricum_genome.fasta.not', 28804),
                                      ('CrossoptilonMantchuricum_genome.fasta.nin', 29000),
                                      ('CrossoptilonMantchuricum_genome.fasta.nto', 9536)],
         'CryptoproctaFerox': [('CryptoproctaFerox_genome.fasta.ntf', 16284),
                               ('CryptoproctaFerox_genome.fasta.nto', 1856048),
                               ('CryptoproctaFerox_genome.fasta.not', 5568340),
                               ('CryptoproctaFerox_genome.fasta.ndb', 20380),
                               ('CryptoproctaFerox_genome.fasta.nsq', 615075730),
                               ('CryptoproctaFerox_genome.fasta.nhr', 79262118),
                               ('CryptoproctaFerox_genome.fasta.nin', 5568520)],
         'FelisCatus': [('FelisCatus_genome.fasta.ntf', 16284),
                        ('FelisCatus_genome.fasta.nto', 21904),
                        ('FelisCatus_genome.fasta.not', 65908),
                        ('FelisCatus_genome.fasta.ndb', 20380),
                        ('FelisCatus_genome.fasta.nsq', 615585257),
                        ('FelisCatus_genome.fasta.nhr', 1106566),
                        ('FelisCatus_genome.fasta.nin', 66080)],
         'GorillaGorilla': [('GorillaGorilla_genome.fasta.ndb', 20380),
                            ('GorillaGorilla_genome.fasta.ntf', 16284),
                            ('GorillaGorilla_genome.fasta.nhr', 1189734),
                            ('GorillaGorilla_genome.fasta.nsq', 761311629),
                            ('GorillaGorilla_genome.fasta.not', 65728),
                            ('GorillaGorilla_genome.fasta.nin', 65908),
                            ('GorillaGorilla_genome.fasta.nto', 21844)],
         'GrammomysDolichurus': [('GrammomysDolichurus_genome.fasta.ndb', 20380),
                                 ('GrammomysDolichurus_genome.fasta.ntf', 16284),
                                 ('GrammomysDolichurus_genome.fasta.nhr', 84331077),
                                 ('GrammomysDolichurus_genome.fasta.nsq', 523427235),
                                 ('GrammomysDolichurus_genome.fasta.not', 6108772),
                                 ('GrammomysDolichurus_genome.fasta.nin', 6108960),
                                 ('GrammomysDolichurus_genome.fasta.nto', 2036192)],
         'GrammomysSurdaster': [('GrammomysSurdaster_genome.fasta.ndb', 20380),
                                ('GrammomysSurdaster_genome.fasta.ntf', 16284),
                                ('GrammomysSurdaster_genome.fasta.nhr', 3566128),
                                ('GrammomysSurdaster_genome.fasta.nsq', 603849422),
                                ('GrammomysSurdaster_genome.fasta.not', 283780),
                                ('GrammomysSurdaster_genome.fasta.nin', 283968),
                                ('GrammomysSurdaster_genome.fasta.nto', 94528)],
         'GuloGulo': [('GuloGulo_genome.fasta.ntf', 16284),
                      ('GuloGulo_genome.fasta.nto', 324180),
                      ('GuloGulo_genome.fasta.not', 972736),
                      ('GuloGulo_genome.fasta.ndb', 20380),
                      ('GuloGulo_genome.fasta.nsq', 563615741),
                      ('GuloGulo_genome.fasta.nhr', 15273820),
                      ('GuloGulo_genome.fasta.nin', 972900)],
         'HyaenaHyaena': [('HyaenaHyaena_genome.fasta.ntf', 16284),
                          ('HyaenaHyaena_genome.fasta.nto', 1400796),
                          ('HyaenaHyaena_genome.fasta.not', 4202584),
                          ('HyaenaHyaena_genome.fasta.ndb', 20380),
                          ('HyaenaHyaena_genome.fasta.nsq', 611724196),
                          ('HyaenaHyaena_genome.fasta.nhr', 58049078),
                          ('HyaenaHyaena_genome.fasta.nin', 4202756)],
         'HylobatesMoloch': [('HylobatesMoloch_genome.fasta.ndb', 20380),
                             ('HylobatesMoloch_genome.fasta.ntf', 16284),
                             ('HylobatesMoloch_genome.fasta.nhr', 2859080),
                             ('HylobatesMoloch_genome.fasta.nsq', 712401066),
                             ('HylobatesMoloch_genome.fasta.not', 220708),
                             ('HylobatesMoloch_genome.fasta.nin', 220888),
                             ('HylobatesMoloch_genome.fasta.nto', 73504)],
         'HylomyscusAlleni': [('HylomyscusAlleni_genome.fasta.ndb', 20380),
                              ('HylomyscusAlleni_genome.fasta.ntf', 16284),
                              ('HylomyscusAlleni_genome.fasta.nhr', 94851497),
                              ('HylomyscusAlleni_genome.fasta.nsq', 546934956),
                              ('HylomyscusAlleni_genome.fasta.not', 7002748),
                              ('HylomyscusAlleni_genome.fasta.nin', 7002928),
                              ('HylomyscusAlleni_genome.fasta.nto', 2334184)],
         'LagopusLeucura': [('LagopusLeucura_genome.fasta.ndb', 20380),
                            ('LagopusLeucura_genome.fasta.ntf', 16284),
                            ('LagopusLeucura_genome.fasta.nhr', 1210309),
                            ('LagopusLeucura_genome.fasta.nsq', 254140926),
                            ('LagopusLeucura_genome.fasta.not', 88024),
                            ('LagopusLeucura_genome.fasta.nin', 88196),
                            ('LagopusLeucura_genome.fasta.nto', 29276)],
         'LagopusMuta': [('LagopusMuta_genome.fasta.ndb', 20380),
                         ('LagopusMuta_genome.fasta.ntf', 16284),
                         ('LagopusMuta_genome.fasta.nin', 2048),
                         ('LagopusMuta_genome.fasta.nhr', 25265),
                         ('LagopusMuta_genome.fasta.nsq', 256690650),
                         ('LagopusMuta_genome.fasta.not', 1876),
                         ('LagopusMuta_genome.fasta.nto', 560)],
         'LophuraNycthemera': [('LophuraNycthemera_genome.fasta.ndb', 20380),
                               ('LophuraNycthemera_genome.fasta.ntf', 16284),
                               ('LophuraNycthemera_genome.fasta.nhr', 356962),
                               ('LophuraNycthemera_genome.fasta.nsq', 253603129),
                               ('LophuraNycthemera_genome.fasta.not', 18544),
                               ('LophuraNycthemera_genome.fasta.nin', 18724),
                               ('LophuraNycthemera_genome.fasta.nto', 6116)],
         'LutraLutra': [('LutraLutra_genome.fasta.ntf', 16284),
                        ('LutraLutra_genome.fasta.nto', 80),
                        ('LutraLutra_genome.fasta.not', 436),
                        ('LutraLutra_genome.fasta.ndb', 20380),
                        ('LutraLutra_genome.fasta.nsq', 609616234),
                        ('LutraLutra_genome.fasta.nhr', 6381),
                        ('LutraLutra_genome.fasta.nin', 608)],
         'LynxRufus': [('LynxRufus_genome.fasta.ntf', 16284),
                       ('LynxRufus_genome.fasta.nto', 208),
                       ('LynxRufus_genome.fasta.not', 820),
                       ('LynxRufus_genome.fasta.ndb', 20380),
                       ('LynxRufus_genome.fasta.nsq', 609814336),
                       ('LynxRufus_genome.fasta.nhr', 11789),
                       ('LynxRufus_genome.fasta.nin', 984)],
         'LyrurusTetrix': [('LyrurusTetrix_genome.fasta.ndb', 20380),
                           ('LyrurusTetrix_genome.fasta.ntf', 16284),
                           ('LyrurusTetrix_genome.fasta.nin', 11664056),
                           ('LyrurusTetrix_genome.fasta.nhr', 144448831),
                           ('LyrurusTetrix_genome.fasta.nsq', 164865265),
                           ('LyrurusTetrix_genome.fasta.not', 11663884),
                           ('LyrurusTetrix_genome.fasta.nto', 3887896)],
         'MacacaMulatta': [('MacacaMulatta_genome.fasta.ndb', 20380),
                           ('MacacaMulatta_genome.fasta.ntf', 16284),
                           ('MacacaMulatta_genome.fasta.nhr', 300185),
                           ('MacacaMulatta_genome.fasta.nsq', 759318537),
                           ('MacacaMulatta_genome.fasta.not', 22480),
                           ('MacacaMulatta_genome.fasta.nin', 22652),
                           ('MacacaMulatta_genome.fasta.nto', 7428)],
         'MastomysCoucha': [('MastomysCoucha_genome.fasta.ndb', 20380),
                            ('MastomysCoucha_genome.fasta.ntf', 16284),
                            ('MastomysCoucha_genome.fasta.nhr', 3962),
                            ('MastomysCoucha_genome.fasta.nsq', 628226218),
                            ('MastomysCoucha_genome.fasta.not', 232),
                            ('MastomysCoucha_genome.fasta.nin', 412),
                            ('MastomysCoucha_genome.fasta.nto', 12)],
         'MastomysNatalensis': [('MastomysNatalensis_genome.fasta.ndb', 20380),
                                ('MastomysNatalensis_genome.fasta.ntf', 16284),
                                ('MastomysNatalensis_genome.fasta.nhr', 52554876),
                                ('MastomysNatalensis_genome.fasta.nsq', 622443336),
                                ('MastomysNatalensis_genome.fasta.not', 3845440),
                                ('MastomysNatalensis_genome.fasta.nin', 3845628),
                                ('MastomysNatalensis_genome.fasta.nto', 1281748)],
         'MeleagrisGallopavo': [('MeleagrisGallopavo_genome.fasta.ndb', 20380),
                                ('MeleagrisGallopavo_genome.fasta.ntf', 16284),
                                ('MeleagrisGallopavo_genome.fasta.nhr', 41380324),
                                ('MeleagrisGallopavo_genome.fasta.nsq', 279551027),
                                ('MeleagrisGallopavo_genome.fasta.not', 2222008),
                                ('MeleagrisGallopavo_genome.fasta.nin', 2222188),
                                ('MeleagrisGallopavo_genome.fasta.nto', 740604)],
         'MelesMeles': [('MelesMeles_genome.fasta.ntf', 16284),
                        ('MelesMeles_genome.fasta.nto', 2056),
                        ('MelesMeles_genome.fasta.not', 6364),
                        ('MelesMeles_genome.fasta.ndb', 20380),
                        ('MelesMeles_genome.fasta.nsq', 684674775),
                        ('MelesMeles_genome.fasta.nhr', 88412),
                        ('MelesMeles_genome.fasta.nin', 6536)],
         'MiroungaLeonina': [('MiroungaLeonina_genome.fasta.ntf', 16284),
                             ('MiroungaLeonina_genome.fasta.nto', 4360),
                             ('MiroungaLeonina_genome.fasta.not', 13276),
                             ('MiroungaLeonina_genome.fasta.ndb', 20380),
                             ('MiroungaLeonina_genome.fasta.nsq', 604534239),
                             ('MiroungaLeonina_genome.fasta.nhr', 174714),
                             ('MiroungaLeonina_genome.fasta.nin', 13456)],
         'MusCaroli': [('MusCaroli_genome.fasta.ndb', 20380),
                       ('MusCaroli_genome.fasta.ntf', 16284),
                       ('MusCaroli_genome.fasta.nhr', 512238),
                       ('MusCaroli_genome.fasta.nsq', 641128168),
                       ('MusCaroli_genome.fasta.not', 37852),
                       ('MusCaroli_genome.fasta.nin', 38016),
                       ('MusCaroli_genome.fasta.nto', 12552)],
         'MusMinutoides': [('MusMinutoides_genome.fasta.ndb', 20380),
                           ('MusMinutoides_genome.fasta.ntf', 16284),
                           ('MusMinutoides_genome.fasta.nhr', 11097558),
                           ('MusMinutoides_genome.fasta.nsq', 709904429),
                           ('MusMinutoides_genome.fasta.not', 796108),
                           ('MusMinutoides_genome.fasta.nin', 796280),
                           ('MusMinutoides_genome.fasta.nto', 265304)],
         'MusPahari': [('MusPahari_genome.fasta.ndb', 20380),
                       ('MusPahari_genome.fasta.ntf', 16284),
                       ('MusPahari_genome.fasta.nhr', 417277),
                       ('MusPahari_genome.fasta.nsq', 621230579),
                       ('MusPahari_genome.fasta.not', 30880),
                       ('MusPahari_genome.fasta.nin', 31044),
                       ('MusPahari_genome.fasta.nto', 10228)],
         'MusSpicilegus': [('MusSpicilegus_genome.fasta.ndb', 20380),
                           ('MusSpicilegus_genome.fasta.ntf', 16284),
                           ('MusSpicilegus_genome.fasta.nhr', 5623836),
                           ('MusSpicilegus_genome.fasta.nsq', 625204940),
                           ('MusSpicilegus_genome.fasta.not', 450340),
                           ('MusSpicilegus_genome.fasta.nin', 450512),
                           ('MusSpicilegus_genome.fasta.nto', 150048)],
         'MusSpretus': [('MusSpretus_genome.fasta.ndb', 20380),
                        ('MusSpretus_genome.fasta.ntf', 16284),
                        ('MusSpretus_genome.fasta.nhr', 827806),
                        ('MusSpretus_genome.fasta.nsq', 659192465),
                        ('MusSpretus_genome.fasta.not', 64756),
                        ('MusSpretus_genome.fasta.nin', 64920),
                        ('MusSpretus_genome.fasta.nto', 21520)],
         'MustelaNigripes': [('MustelaNigripes_genome.fasta.ntf', 16284),
                             ('MustelaNigripes_genome.fasta.nto', 82200),
                             ('MustelaNigripes_genome.fasta.not', 246796),
                             ('MustelaNigripes_genome.fasta.ndb', 20380),
                             ('MustelaNigripes_genome.fasta.nsq', 625051597),
                             ('MustelaNigripes_genome.fasta.nhr', 3342077),
                             ('MustelaNigripes_genome.fasta.nin', 246976)],
         'NomascusLeucogenys': [('NomascusLeucogenys_genome.fasta.ndb', 20380),
                                ('NomascusLeucogenys_genome.fasta.ntf', 16284),
                                ('NomascusLeucogenys_genome.fasta.nhr', 360108),
                                ('NomascusLeucogenys_genome.fasta.nsq', 711025689),
                                ('NomascusLeucogenys_genome.fasta.not', 26920),
                                ('NomascusLeucogenys_genome.fasta.nin', 27108),
                                ('NomascusLeucogenys_genome.fasta.nto', 8908)],
         'OdobenusRosmarus': [('OdobenusRosmarus_genome.fasta.ntf', 16284),
                              ('OdobenusRosmarus_genome.fasta.nto', 15472),
                              ('OdobenusRosmarus_genome.fasta.not', 46612),
                              ('OdobenusRosmarus_genome.fasta.ndb', 20380),
                              ('OdobenusRosmarus_genome.fasta.nsq', 600644137),
                              ('OdobenusRosmarus_genome.fasta.nhr', 692615),
                              ('OdobenusRosmarus_genome.fasta.nin', 46792)],
         'PanTroglodytes': [('PanTroglodytes_genome.fasta.ndb', 20380),
                            ('PanTroglodytes_genome.fasta.ntf', 16284),
                            ('PanTroglodytes_genome.fasta.nhr', 868195),
                            ('PanTroglodytes_genome.fasta.nsq', 756074970),
                            ('PanTroglodytes_genome.fasta.not', 52036),
                            ('PanTroglodytes_genome.fasta.nin', 52216),
                            ('PanTroglodytes_genome.fasta.nto', 17280)],
         'PantheraTigris': [('PantheraTigris_genome.fasta.nhr', 200683),
                            ('PantheraTigris_genome.fasta.ntf', 16284),
                            ('PantheraTigris_genome.fasta.nto', 4712),
                            ('PantheraTigris_genome.fasta.not', 14332),
                            ('PantheraTigris_genome.fasta.ndb', 20380),
                            ('PantheraTigris_genome.fasta.nsq', 602112451),
                            ('PantheraTigris_genome.fasta.nin', 14512)],
         'PapioAnubis': [('PapioAnubis_genome.fasta.ndb', 20380),
                         ('PapioAnubis_genome.fasta.ntf', 16284),
                         ('PapioAnubis_genome.fasta.nhr', 1660831),
                         ('PapioAnubis_genome.fasta.nsq', 717494436),
                         ('PapioAnubis_genome.fasta.not', 133636),
                         ('PapioAnubis_genome.fasta.nin', 133808),
                         ('PapioAnubis_genome.fasta.nto', 44480)],
         'ParadoxurusHermaphroditus': [('ParadoxurusHermaphroditus_genome.fasta.ntf',
                                        16284),
                                       ('ParadoxurusHermaphroditus_genome.fasta.nto', 1752404),
                                       ('ParadoxurusHermaphroditus_genome.fasta.not', 5257408),
                                       ('ParadoxurusHermaphroditus_genome.fasta.ndb', 20380),
                                       ('ParadoxurusHermaphroditus_genome.fasta.nsq', 619630145),
                                       ('ParadoxurusHermaphroditus_genome.fasta.nhr', 78774558),
                                       ('ParadoxurusHermaphroditus_genome.fasta.nin', 5257604)],
         'PavoCristatus': [('PavoCristatus_genome.fasta.ndb', 20380),
                           ('PavoCristatus_genome.fasta.ntf', 16284),
                           ('PavoCristatus_genome.fasta.nhr', 110027),
                           ('PavoCristatus_genome.fasta.nsq', 261688314),
                           ('PavoCristatus_genome.fasta.not', 8620),
                           ('PavoCristatus_genome.fasta.nin', 8792),
                           ('PavoCristatus_genome.fasta.nto', 2808)],
         'PavoMuticus': [('PavoMuticus_genome.fasta.ndb', 20380),
                         ('PavoMuticus_genome.fasta.ntf', 16284),
                         ('PavoMuticus_genome.fasta.nin', 29432),
                         ('PavoMuticus_genome.fasta.nhr', 375349),
                         ('PavoMuticus_genome.fasta.nsq', 265553754),
                         ('PavoMuticus_genome.fasta.not', 29260),
                         ('PavoMuticus_genome.fasta.nto', 9688)],
         'PhasianusColchicus': [('PhasianusColchicus_genome.fasta.ndb', 20380),
                                ('PhasianusColchicus_genome.fasta.ntf', 16284),
                                ('PhasianusColchicus_genome.fasta.nhr', 6665760),
                                ('PhasianusColchicus_genome.fasta.nsq', 254584785),
                                ('PhasianusColchicus_genome.fasta.not', 476032),
                                ('PhasianusColchicus_genome.fasta.nin', 476212),
                                ('PhasianusColchicus_genome.fasta.nto', 158612)],
         'PiliocolobusTephrosceles': [('PiliocolobusTephrosceles_genome.fasta.ndb',
                                       20380),
                                      ('PiliocolobusTephrosceles_genome.fasta.ntf', 16284),
                                      ('PiliocolobusTephrosceles_genome.fasta.nhr', 7840206),
                                      ('PiliocolobusTephrosceles_genome.fasta.nsq', 760211459),
                                      ('PiliocolobusTephrosceles_genome.fasta.not', 559756),
                                      ('PiliocolobusTephrosceles_genome.fasta.nin', 559952),
                                      ('PiliocolobusTephrosceles_genome.fasta.nto', 186520)],
         'PitheciaPithecia': [('PitheciaPithecia_genome.fasta.ndb', 20380),
                              ('PitheciaPithecia_genome.fasta.ntf', 16284),
                              ('PitheciaPithecia_genome.fasta.nhr', 146647250),
                              ('PitheciaPithecia_genome.fasta.nsq', 734089273),
                              ('PitheciaPithecia_genome.fasta.not', 10707868),
                              ('PitheciaPithecia_genome.fasta.nin', 10708048),
                              ('PitheciaPithecia_genome.fasta.nto', 3569224)],
         'PlecturocebusDonacophilus': [('PlecturocebusDonacophilus_genome.fasta.ndb',
                                        20380),
                                       ('PlecturocebusDonacophilus_genome.fasta.ntf', 16284),
                                       ('PlecturocebusDonacophilus_genome.fasta.nhr', 180360420),
                                       ('PlecturocebusDonacophilus_genome.fasta.nsq', 740680183),
                                       ('PlecturocebusDonacophilus_genome.fasta.not', 12407704),
                                       ('PlecturocebusDonacophilus_genome.fasta.nin', 12407900),
                                       ('PlecturocebusDonacophilus_genome.fasta.nto', 4135836)],
         'PongoAbelli': [('PongoAbelli_genome.fasta.ndb', 20380),
                         ('PongoAbelli_genome.fasta.ntf', 16284),
                         ('PongoAbelli_genome.fasta.nhr', 866693),
                         ('PongoAbelli_genome.fasta.nsq', 766308041),
                         ('PongoAbelli_genome.fasta.not', 63028),
                         ('PongoAbelli_genome.fasta.nin', 63200),
                         ('PongoAbelli_genome.fasta.nto', 20944)],
         'PraomysDelectorum': [('PraomysDelectorum_genome.fasta.ndb', 20380),
                               ('PraomysDelectorum_genome.fasta.ntf', 16284),
                               ('PraomysDelectorum_genome.fasta.nhr', 40269687),
                               ('PraomysDelectorum_genome.fasta.nsq', 602535043),
                               ('PraomysDelectorum_genome.fasta.not', 2967592),
                               ('PraomysDelectorum_genome.fasta.nin', 2967772),
                               ('PraomysDelectorum_genome.fasta.nto', 989132)],
         'ProcyonLotor': [('ProcyonLotor_genome.fasta.ntf', 16284),
                          ('ProcyonLotor_genome.fasta.nto', 196908),
                          ('ProcyonLotor_genome.fasta.not', 590920),
                          ('ProcyonLotor_genome.fasta.ndb', 20380),
                          ('ProcyonLotor_genome.fasta.nsq', 565200128),
                          ('ProcyonLotor_genome.fasta.nhr', 6999179),
                          ('ProcyonLotor_genome.fasta.nin', 591092)],
         'RattusNorvegicus': [('RattusNorvegicus_genome.fasta.ndb', 20380),
                              ('RattusNorvegicus_genome.fasta.ntf', 16284),
                              ('RattusNorvegicus_genome.fasta.nhr', 164415),
                              ('RattusNorvegicus_genome.fasta.nsq', 718887184),
                              ('RattusNorvegicus_genome.fasta.not', 11368),
                              ('RattusNorvegicus_genome.fasta.nin', 11548),
                              ('RattusNorvegicus_genome.fasta.nto', 3724)],
         'RattusRattus': [('RattusRattus_genome.fasta.ndb', 20380),
                          ('RattusRattus_genome.fasta.ntf', 16284),
                          ('RattusRattus_genome.fasta.nhr', 363442),
                          ('RattusRattus_genome.fasta.nsq', 595743366),
                          ('RattusRattus_genome.fasta.not', 25972),
                          ('RattusRattus_genome.fasta.nin', 26144),
                          ('RattusRattus_genome.fasta.nto', 8592)],
         'RhabdomysDilectus': [('RhabdomysDilectus_genome.fasta.ndb', 20380),
                               ('RhabdomysDilectus_genome.fasta.ntf', 16284),
                               ('RhabdomysDilectus_genome.fasta.nhr', 4877327),
                               ('RhabdomysDilectus_genome.fasta.nsq', 565469982),
                               ('RhabdomysDilectus_genome.fasta.not', 371344),
                               ('RhabdomysDilectus_genome.fasta.nin', 371524),
                               ('RhabdomysDilectus_genome.fasta.nto', 123716)],
         'RhynchomysSoricoides': [('RhynchomysSoricoides_genome.fasta.ndb', 20380),
                                  ('RhynchomysSoricoides_genome.fasta.ntf', 16284),
                                  ('RhynchomysSoricoides_genome.fasta.nhr', 2273649),
                                  ('RhynchomysSoricoides_genome.fasta.nsq', 546256123),
                                  ('RhynchomysSoricoides_genome.fasta.not', 170224),
                                  ('RhynchomysSoricoides_genome.fasta.nin', 170412),
                                  ('RhynchomysSoricoides_genome.fasta.nto', 56676)],
         'SaimiriBoliviensis': [('SaimiriBoliviensis_genome.fasta.ndb', 20380),
                                ('SaimiriBoliviensis_genome.fasta.ntf', 16284),
                                ('SaimiriBoliviensis_genome.fasta.nhr', 349608),
                                ('SaimiriBoliviensis_genome.fasta.nsq', 662877877),
                                ('SaimiriBoliviensis_genome.fasta.not', 20488),
                                ('SaimiriBoliviensis_genome.fasta.nin', 20668),
                                ('SaimiriBoliviensis_genome.fasta.nto', 6764)],
         'SpeothosVenaticus': [('SpeothosVenaticus_genome.fasta.ntf', 16284),
                               ('SpeothosVenaticus_genome.fasta.nto', 210316),
                               ('SpeothosVenaticus_genome.fasta.not', 631144),
                               ('SpeothosVenaticus_genome.fasta.ndb', 20380),
                               ('SpeothosVenaticus_genome.fasta.nsq', 580856849),
                               ('SpeothosVenaticus_genome.fasta.nhr', 7846348),
                               ('SpeothosVenaticus_genome.fasta.nin', 631324)],
         'SpilogaleGracilis': [('SpilogaleGracilis_genome.fasta.ntf', 16284),
                               ('SpilogaleGracilis_genome.fasta.nto', 1687940),
                               ('SpilogaleGracilis_genome.fasta.not', 5064016),
                               ('SpilogaleGracilis_genome.fasta.ndb', 20380),
                               ('SpilogaleGracilis_genome.fasta.nsq', 616229541),
                               ('SpilogaleGracilis_genome.fasta.nhr', 72497726),
                               ('SpilogaleGracilis_genome.fasta.nin', 5064196)],
         'SuricataSuricatta': [('SuricataSuricatta_genome.fasta.ntf', 16284),
                               ('SuricataSuricatta_genome.fasta.nto', 255068),
                               ('SuricataSuricatta_genome.fasta.not', 765400),
                               ('SuricataSuricatta_genome.fasta.ndb', 20380),
                               ('SuricataSuricatta_genome.fasta.nsq', 588957113),
                               ('SuricataSuricatta_genome.fasta.nhr', 10417513),
                               ('SuricataSuricatta_genome.fasta.nin', 765580)],
         'SyrmaticusMikado': [('SyrmaticusMikado_genome.fasta.ndb', 20380),
                              ('SyrmaticusMikado_genome.fasta.ntf', 16284),
                              ('SyrmaticusMikado_genome.fasta.nhr', 11147017),
                              ('SyrmaticusMikado_genome.fasta.nsq', 267358598),
                              ('SyrmaticusMikado_genome.fasta.not', 855280),
                              ('SyrmaticusMikado_genome.fasta.nin', 855460),
                              ('SyrmaticusMikado_genome.fasta.nto', 285028)],
         'TrachypithecusFrancoisi': [('TrachypithecusFrancoisi_genome.fasta.ndb',
                                      20380),
                                     ('TrachypithecusFrancoisi_genome.fasta.ntf', 16284),
                                     ('TrachypithecusFrancoisi_genome.fasta.nhr', 811475),
                                     ('TrachypithecusFrancoisi_genome.fasta.nsq', 725919335),
                                     ('TrachypithecusFrancoisi_genome.fasta.not', 53272),
                                     ('TrachypithecusFrancoisi_genome.fasta.nin', 53468),
                                     ('TrachypithecusFrancoisi_genome.fasta.nto', 17692)],
         'TympanuchusCupido': [('TympanuchusCupido_genome.fasta.ndb', 20380),
                               ('TympanuchusCupido_genome.fasta.ntf', 16284),
                               ('TympanuchusCupido_genome.fasta.nhr', 2056961),
                               ('TympanuchusCupido_genome.fasta.nsq', 246353243),
                               ('TympanuchusCupido_genome.fasta.not', 146140),
                               ('TympanuchusCupido_genome.fasta.nin', 146320),
                               ('TympanuchusCupido_genome.fasta.nto', 48648)],
         'UrsusAmericanus': [('UrsusAmericanus_genome.fasta.ntf', 16284),
                             ('UrsusAmericanus_genome.fasta.nto', 445884),
                             ('UrsusAmericanus_genome.fasta.not', 1337848),
                             ('UrsusAmericanus_genome.fasta.ndb', 20380),
                             ('UrsusAmericanus_genome.fasta.nsq', 650150246),
                             ('UrsusAmericanus_genome.fasta.nhr', 17585509),
                             ('UrsusAmericanus_genome.fasta.nin', 1338028)],
         'UrsusMaritimus': [('UrsusMaritimus_genome.fasta.ntf', 16284),
                            ('UrsusMaritimus_genome.fasta.nto', 15504),
                            ('UrsusMaritimus_genome.fasta.not', 46708),
                            ('UrsusMaritimus_genome.fasta.ndb', 20380),
                            ('UrsusMaritimus_genome.fasta.nsq', 582998576),
                            ('UrsusMaritimus_genome.fasta.nhr', 603451),
                            ('UrsusMaritimus_genome.fasta.nin', 46888)],
         'VulpesFerrilata': [('VulpesFerrilata_genome.fasta.ntf', 16284),
                             ('VulpesFerrilata_genome.fasta.nto', 840),
                             ('VulpesFerrilata_genome.fasta.not', 2716),
                             ('VulpesFerrilata_genome.fasta.ndb', 20380),
                             ('VulpesFerrilata_genome.fasta.nsq', 594915962),
                             ('VulpesFerrilata_genome.fasta.nhr', 35871),
                             ('VulpesFerrilata_genome.fasta.nin', 2896)],
         'ZalophusCalifornianus': [('ZalophusCalifornianus_genome.fasta.ntf', 16284),
                                   ('ZalophusCalifornianus_genome.fasta.nto', 92),
                                   ('ZalophusCalifornianus_genome.fasta.not', 472),
                                   ('ZalophusCalifornianus_genome.fasta.ndb', 20380),
                                   ('ZalophusCalifornianus_genome.fasta.nsq', 602458774),
                                   ('ZalophusCalifornianus_genome.fasta.nhr', 7529),
                                   ('ZalophusCalifornianus_genome.fasta.nin', 660)],
         'DrosophilaSimulans': [('DrosophilaSimulans_genome.fasta.ndb', 20380),
                                ('DrosophilaSimulans_genome.fasta.ntf', 16284),
                                ('DrosophilaSimulans_genome.fasta.nhr', 41496),
                                ('DrosophilaSimulans_genome.fasta.nsq', 33392323),
                                ('DrosophilaSimulans_genome.fasta.not', 2932),
                                ('DrosophilaSimulans_genome.fasta.nin', 3112),
                                ('DrosophilaSimulans_genome.fasta.nto', 912)],
         'DrosophilaMauritiana': [('DrosophilaMauritiana_genome.fasta.ndb', 20380),
                                  ('DrosophilaMauritiana_genome.fasta.ntf', 16284),
                                  ('DrosophilaMauritiana_genome.fasta.nhr', 78989),
                                  ('DrosophilaMauritiana_genome.fasta.nsq', 32533720),
                                  ('DrosophilaMauritiana_genome.fasta.not', 5608),
                                  ('DrosophilaMauritiana_genome.fasta.nin', 5788),
                                  ('DrosophilaMauritiana_genome.fasta.nto', 1804)],
         'DrosophilaSechellia': [('DrosophilaSechellia_genome.fasta.ndb', 20380),
                                 ('DrosophilaSechellia_genome.fasta.ntf', 16284),
                                 ('DrosophilaSechellia_genome.fasta.nhr', 88806),
                                 ('DrosophilaSechellia_genome.fasta.nsq', 36976228),
                                 ('DrosophilaSechellia_genome.fasta.not', 6352),
                                 ('DrosophilaSechellia_genome.fasta.nin', 6532),
                                 ('DrosophilaSechellia_genome.fasta.nto', 2052)],
         'DrosophilaYakuba': [('DrosophilaYakuba_genome.fasta.ndb', 20380),
                              ('DrosophilaYakuba_genome.fasta.ntf', 16284),
                              ('DrosophilaYakuba_genome.fasta.nhr', 62908),
                              ('DrosophilaYakuba_genome.fasta.nsq', 37827741),
                              ('DrosophilaYakuba_genome.fasta.not', 4564),
                              ('DrosophilaYakuba_genome.fasta.nin', 4736),
                              ('DrosophilaYakuba_genome.fasta.nto', 1456)],
         'DrosophilaSantomea': [('DrosophilaSantomea_genome.fasta.ndb', 20380),
                               ('DrosophilaSantomea_genome.fasta.ntf', 16284),
                               ('DrosophilaSantomea_genome.fasta.nhr', 18560),
                               ('DrosophilaSantomea_genome.fasta.nsq', 36704225),
                               ('DrosophilaSantomea_genome.fasta.not', 1168),
                               ('DrosophilaSantomea_genome.fasta.nin', 1348),
                               ('DrosophilaSantomea_genome.fasta.nto', 324)],
         'DrosophilaTeissieri': [('DrosophilaTeissieri_genome.fasta.ndb', 20380),
                                 ('DrosophilaTeissieri_genome.fasta.ntf', 16284),
                                 ('DrosophilaTeissieri_genome.fasta.nhr', 179447),
                                 ('DrosophilaTeissieri_genome.fasta.nsq', 36173752),
                                 ('DrosophilaTeissieri_genome.fasta.not', 13636),
                                 ('DrosophilaTeissieri_genome.fasta.nin', 13816),
                                 ('DrosophilaTeissieri_genome.fasta.nto', 4480)],
         'DrosophilaOrena': [('DrosophilaOrena_genome.fasta.ndb', 20380),
                            ('DrosophilaOrena_genome.fasta.ntf', 16284),
                            ('DrosophilaOrena_genome.fasta.nhr', 74388),
                            ('DrosophilaOrena_genome.fasta.nsq', 45723031),
                            ('DrosophilaOrena_genome.fasta.not', 4972),
                            ('DrosophilaOrena_genome.fasta.nin', 5144),
                            ('DrosophilaOrena_genome.fasta.nto', 1592)],
         'DrosophilaErecta': [('DrosophilaErecta_genome.fasta.ndb', 20380),
                              ('DrosophilaErecta_genome.fasta.ntf', 16284),
                              ('DrosophilaErecta_genome.fasta.nhr', 81003),
                              ('DrosophilaErecta_genome.fasta.nsq', 33914625),
                              ('DrosophilaErecta_genome.fasta.not', 5896),
                              ('DrosophilaErecta_genome.fasta.nin', 6068),
                              ('DrosophilaErecta_genome.fasta.nto', 1900)],
         'DrosophilaEugracilis': [('DrosophilaEugracilis_genome.fasta.ndb', 20380),
                                  ('DrosophilaEugracilis_genome.fasta.ntf', 16284),
                                  ('DrosophilaEugracilis_genome.fasta.nhr', 361216),
                                  ('DrosophilaEugracilis_genome.fasta.nsq', 41204062),
                                  ('DrosophilaEugracilis_genome.fasta.not', 25804),
                                  ('DrosophilaEugracilis_genome.fasta.nin', 25984),
                                  ('DrosophilaEugracilis_genome.fasta.nto', 8536)],
         'DrosophilaSubpulchrella': [('DrosophilaSubpulchrella_genome.fasta.ndb',
                                      20380),
                                     ('DrosophilaSubpulchrella_genome.fasta.ntf', 16284),
                                     ('DrosophilaSubpulchrella_genome.fasta.nhr', 220613),
                                     ('DrosophilaSubpulchrella_genome.fasta.nsq', 67022015),
                                     ('DrosophilaSubpulchrella_genome.fasta.not', 16552),
                                     ('DrosophilaSubpulchrella_genome.fasta.nin', 16740),
                                     ('DrosophilaSubpulchrella_genome.fasta.nto', 5452)],
         'DrosophilaSuzukii': [('DrosophilaSuzukii_genome.fasta.ndb', 20380),
                                    ('DrosophilaSuzukii_genome.fasta.ntf', 16284),
                                    ('DrosophilaSuzukii_genome.fasta.nhr', 84435),
                                    ('DrosophilaSuzukii_genome.fasta.nsq', 67003287),
                                    ('DrosophilaSuzukii_genome.fasta.not', 6460),
                                    ('DrosophilaSuzukii_genome.fasta.nin', 6640),
                                    ('DrosophilaSuzukii_genome.fasta.nto', 2088)],
         'DrosophilaBiarmipes': [('DrosophilaBiarmipes_genome.fasta.ndb', 20380),
                                 ('DrosophilaBiarmipes_genome.fasta.ntf', 16284),
                                 ('DrosophilaBiarmipes_genome.fasta.nhr', 46642),
                                 ('DrosophilaBiarmipes_genome.fasta.nsq', 46329699),
                                 ('DrosophilaBiarmipes_genome.fasta.not', 3304),
                                 ('DrosophilaBiarmipes_genome.fasta.nin', 3484),
                                 ('DrosophilaBiarmipes_genome.fasta.nto', 1036)],
         'DrosophilaTakahashii': [('DrosophilaTakahashii_genome.fasta.ndb', 20380),
                                  ('DrosophilaTakahashii_genome.fasta.ntf', 16284),
                                  ('DrosophilaTakahashii_genome.fasta.nhr', 20166),
                                  ('DrosophilaTakahashii_genome.fasta.nsq', 41382067),
                                  ('DrosophilaTakahashii_genome.fasta.not', 1372),
                                  ('DrosophilaTakahashii_genome.fasta.nin', 1552),
                                  ('DrosophilaTakahashii_genome.fasta.nto', 392)],
         'DrosophilaFicusphila': [('DrosophilaFicusphila_genome.fasta.ndb', 20380),
                                  ('DrosophilaFicusphila_genome.fasta.ntf', 16284),
                                  ('DrosophilaFicusphila_genome.fasta.nhr', 139610),
                                  ('DrosophilaFicusphila_genome.fasta.nsq', 41958669),
                                  ('DrosophilaFicusphila_genome.fasta.not', 9964),
                                  ('DrosophilaFicusphila_genome.fasta.nin', 10144),
                                  ('DrosophilaFicusphila_genome.fasta.nto', 3256)],
         'DrosophilaCarrolli': [('DrosophilaCarrolli_genome.fasta.ndb', 20380),
                                ('DrosophilaCarrolli_genome.fasta.ntf', 16284),
                                ('DrosophilaCarrolli_genome.fasta.nhr', 52730),
                                ('DrosophilaCarrolli_genome.fasta.nsq', 57804927),
                                ('DrosophilaCarrolli_genome.fasta.not', 3964),
                                ('DrosophilaCarrolli_genome.fasta.nin', 4144),
                                ('DrosophilaCarrolli_genome.fasta.nto', 1256)],
         'DrosophilaRhopaloa': [('DrosophilaRhopaloa_genome.fasta.ndb', 20380),
                                ('DrosophilaRhopaloa_genome.fasta.ntf', 16284),
                                ('DrosophilaRhopaloa_genome.fasta.nhr', 37284),
                                ('DrosophilaRhopaloa_genome.fasta.nsq', 48377100),
                                ('DrosophilaRhopaloa_genome.fasta.not', 2644),
                                ('DrosophilaRhopaloa_genome.fasta.nin', 2824),
                                ('DrosophilaRhopaloa_genome.fasta.nto', 816)],
         'DrosophilaKurseongensis': [('DrosophilaKurseongensis_genome.fasta.ndb',
                                      20380),
                                     ('DrosophilaKurseongensis_genome.fasta.ntf', 16284),
                                     ('DrosophilaKurseongensis_genome.fasta.nhr', 43185),
                                     ('DrosophilaKurseongensis_genome.fasta.nsq', 51610268),
                                     ('DrosophilaKurseongensis_genome.fasta.not', 3112),
                                     ('DrosophilaKurseongensis_genome.fasta.nin', 3300),
                                     ('DrosophilaKurseongensis_genome.fasta.nto', 972)],
         'DrosophilaFuyamai': [('DrosophilaFuyamai_genome.fasta.ndb', 20380),
                               ('DrosophilaFuyamai_genome.fasta.ntf', 16284),
                               ('DrosophilaFuyamai_genome.fasta.nhr', 145460),
                               ('DrosophilaFuyamai_genome.fasta.nsq', 57275609),
                               ('DrosophilaFuyamai_genome.fasta.not', 10576),
                               ('DrosophilaFuyamai_genome.fasta.nin', 10748),
                               ('DrosophilaFuyamai_genome.fasta.nto', 3460)],
         'DrosophilaElegans': [('DrosophilaElegans_genome.fasta.ndb', 20380),
                               ('DrosophilaElegans_genome.fasta.ntf', 16284),
                               ('DrosophilaElegans_genome.fasta.nhr', 117792),
                               ('DrosophilaElegans_genome.fasta.nsq', 44611641),
                               ('DrosophilaElegans_genome.fasta.not', 8548),
                               ('DrosophilaElegans_genome.fasta.nin', 8720),
                               ('DrosophilaElegans_genome.fasta.nto', 2784)],
         'PanPaniscus': [('PanPaniscus_genome.fasta.ndb', 20380),
                          ('PanPaniscus_genome.fasta.ntf', 16284),
                          ('PanPaniscus_genome.fasta.nhr', 63393),
                          ('PanPaniscus_genome.fasta.nsq', 817912060),
                          ('PanPaniscus_genome.fasta.not', 5236),
                          ('PanPaniscus_genome.fasta.nto', 1680),
                          ('PanPaniscus_genome.fasta.nin', 5400)],
        'PongoPygmaeus': [('PongoPygmaeus_genome.fasta.ndb', 20380),
                          ('PongoPygmaeus_genome.fasta.ntf', 16284),
                          ('PongoPygmaeus_genome.fasta.nhr', 251405),
                          ('PongoPygmaeus_genome.fasta.nsq', 812347770),
                          ('PongoPygmaeus_genome.fasta.not', 18568),
                          ('PongoPygmaeus_genome.fasta.nin', 18732),
                          ('PongoPygmaeus_genome.fasta.nto', 6124)],
         'NomascusSiki': [('NomascusSiki_genome.fasta.ndb', 20380),
                          ('NomascusSiki_genome.fasta.ntf', 16284),
                          ('NomascusSiki_genome.fasta.nhr', 416888),
                          ('NomascusSiki_genome.fasta.nsq', 695785408),
                          ('NomascusSiki_genome.fasta.not', 31672),
                          ('NomascusSiki_genome.fasta.nin', 31836),
                          ('NomascusSiki_genome.fasta.nto', 10492)],
         'HylobatesPileatus': [('HylobatesPileatus_genome.fasta.ndb', 20380),
                                  ('HylobatesPileatus_genome.fasta.ntf', 16284),
                                  ('HylobatesPileatus_genome.fasta.nhr', 125820),
                                  ('HylobatesPileatus_genome.fasta.nsq', 712857374),
                                  ('HylobatesPileatus_genome.fasta.not', 8824),
                                  ('HylobatesPileatus_genome.fasta.nin', 8996),
                                  ('HylobatesPileatus_genome.fasta.nto', 2876)],
         'SymphalangusSyndactylus': [('SymphalangusSyndactylus_genome.fasta.ndb', 20380),
                                      ('SymphalangusSyndactylus_genome.fasta.ntf', 16284),
                                      ('SymphalangusSyndactylus_genome.fasta.nhr', 102476),
                                      ('SymphalangusSyndactylus_genome.fasta.nsq', 795735238),
                                      ('SymphalangusSyndactylus_genome.fasta.not', 7180),
                                      ('SymphalangusSyndactylus_genome.fasta.nin', 7368),
                                      ('SymphalangusSyndactylus_genome.fasta.nto', 2328)],
         'HoolockLeuconedys': [('HoolockLeuconedys_genome.fasta.ndb', 20380),
                                  ('HoolockLeuconedys_genome.fasta.ntf', 16284),
                                  ('HoolockLeuconedys_genome.fasta.nhr', 467768),
                                  ('HoolockLeuconedys_genome.fasta.nsq', 696251645),
                                  ('HoolockLeuconedys_genome.fasta.not', 35536),
                                  ('HoolockLeuconedys_genome.fasta.nin', 35708),
                                  ('HoolockLeuconedys_genome.fasta.nto', 11780)],
         'OtocolobusManul': [('OtocolobusManul_genome.fasta.ndb', 0),
                                     ('OtocolobusManul_genome.fasta.ntf', 0),
                                     ('OtocolobusManul_genome.fasta.nhr', 0),
                                     ('OtocolobusManul_genome.fasta.nsq', 0),
                                     ('OtocolobusManul_genome.fasta.not', 0),
                                     ('OtocolobusManul_genome.fasta.nin', 0),
                                     ('OtocolobusManul_genome.fasta.nto', 0)],
         'PrionailurusViverrinus': [('PrionailurusViverrinus_genome.fasta.ndb', 0),
                                   ('PrionailurusViverrinus_genome.fasta.ntf', 0),
                                   ('PrionailurusViverrinus_genome.fasta.nhr', 0),
                                   ('PrionailurusViverrinus_genome.fasta.nsq', 0),
                                   ('PrionailurusViverrinus_genome.fasta.not', 0),
                                   ('PrionailurusViverrinus_genome.fasta.nin', 0),
                                   ('PrionailurusViverrinus_genome.fasta.nto', 0)],
         'AcinonyxJubatus': [('AcinonyxJubatus_genome.fasta.ndb', 0),
                                          ('AcinonyxJubatus_genome.fasta.ntf', 0),
                                          ('AcinonyxJubatus_genome.fasta.nhr', 0),
                                          ('AcinonyxJubatus_genome.fasta.nsq', 0),
                                          ('AcinonyxJubatus_genome.fasta.not', 0),
                                          ('AcinonyxJubatus_genome.fasta.nin', 0),
                                          ('AcinonyxJubatus_genome.fasta.nto', 0)],
         'PumaConcolor': [('PumaConcolor_genome.fasta.ndb', 0),
                                   ('PumaConcolor_genome.fasta.ntf', 0),
                                   ('PumaConcolor_genome.fasta.nhr', 0),
                                   ('PumaConcolor_genome.fasta.nsq', 0),
                                   ('PumaConcolor_genome.fasta.not', 0),
                                   ('PumaConcolor_genome.fasta.nin', 0),
                                   ('PumaConcolor_genome.fasta.nto', 0)],
         'NeofelisNebulosa': [('NeofelisNebulosa_genome.fasta.ndb', 0),
                                ('NeofelisNebulosa_genome.fasta.ntf', 0),
                                ('NeofelisNebulosa_genome.fasta.nhr', 0),
                                ('NeofelisNebulosa_genome.fasta.nsq', 0),
                                ('NeofelisNebulosa_genome.fasta.not', 0),
                                ('NeofelisNebulosa_genome.fasta.nin', 0),
                                ('NeofelisNebulosa_genome.fasta.nto', 0)],
         'MungosMungo': [('MungosMungo_genome.fasta.ndb', 0),
                                    ('MungosMungo_genome.fasta.ntf', 0),
                                    ('MungosMungo_genome.fasta.nhr', 0),
                                    ('MungosMungo_genome.fasta.nsq', 0),
                                    ('MungosMungo_genome.fasta.not', 0),
                                    ('MungosMungo_genome.fasta.nin', 0),
                                    ('MungosMungo_genome.fasta.nto', 0)]}

    # define the expected sizes of all possible databases
    expected_databases = {}
    for file in possible_databases[genome.replace("_genome", "")]:
        expected_databases[file[0]] = file[1]

    # check actual size of each expected database
    generated_databases = {}
    for subdirs, dirs, files in os.walk(database_path):
        for file in files:
            # make sure to avoid hidden files
            if file in expected_databases and not file.startswith("."):
                generated_databases[file] = os.path.getsize(database_path + file)

    # no databases have been detected
    if len(generated_databases) == 0:
        message = "\nBlast databases are missing for: %s" % genome
        logging.info(message)
        return False

    # check if all databases are present
    if len(expected_databases) != len(generated_databases):
        message = "\n...WARNING... : Genome : %s blast databases failed to build (likely interrupted)\n" % genome
        logging.info(message)

        # remove these databases
        for file in generated_databases:
            os.remove(database_path + file)
        return False

    # check sizes
    for file in generated_databases:
        if file in expected_databases:
            if os.path.getsize(database_path + file) < expected_databases[file]:
                message = "\n...WARNING... : Genome : %s blast databases were built partially (likely interrupted)\n" \
                                        % genome
                logging.info(message)

                # remove these databases
                for f in generated_databases:
                    os.remove(database_path + f)
                return False

    # if all expected databases are accounted for and of correct sizes
    message = "\nGenome : %s blast database is present" % genome
    logging.info(message)
    return True














