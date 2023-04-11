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

Indexes a genome of interest using a BioPython module SeqIO

"""

from Bio import SeqIO
import logging
import glob
import os


def index_genome_database(wdir, genome_name):
    """Indexes genome databases"""
    
    genomes_dir = wdir + "Genomes/"
    
    name = glob.glob(genomes_dir + genome_name + ".fasta")
    index_path = glob.glob(genomes_dir + genome_name + ".idx")

    if not index_path:
        # it will make a new index
        message = "\nIndexing genome: " + genome_name + "\n"
        print(message)
        logging.info(message)
        genome_index = SeqIO.index_db(genomes_dir + genome_name + ".idx", name[0], "fasta")
        message = "\nGenome: " + genome_name + " has been indexed.\n"
        print(message)
        logging.info(message)

    if index_path:
        # check if index is not partial
        index_complete = check_index(index_path[0])

        if index_complete:
            # it will find the index and put it in the genome_index variable
            genome_index = SeqIO.index_db(genomes_dir + genome_name + ".idx", name[0], "fasta")
            message = "\nGenome: " + genome_name + " index already exists.\n"
            print(message)
            logging.info(message)

        else:
            message = "\n...WARNING... : Genome: " + genome_name + " has been only partially indexed. \
                        \n     -> removing the index and indexing again..."
            print(message)
            logging.info(message)
            os.remove(index_path[0])
            # it will make a new index
            genome_index = SeqIO.index_db(genomes_dir + genome_name + ".idx", name[0], "fasta")
            message = "\nGenome: " + genome_name + " has been indexed.\n"
            print(message)
            logging.info(message)
        
    return genome_index


def check_index(index_path):
    """Checks if indexing was successful and there are no partial indexes"""

    # -5kb for each index
    possible_indexes = {
         'AilurusFulgens_genome.idx': [655355],
         'AlectorisRufa_genome.idx': [643067],
         'AlouattaPalliata_genome.idx': [66265083],
         'AotusNancymaae_genome.idx': [1458171],
         'ApodemusSpeciosus_genome.idx': [18956283],
         'ApodemusSylvaticus_genome.idx': [31600635],
         'ArvicanthisNiloticus_genome.idx': [118779],
         'AtelesGeoffroyi_genome.idx': [49098747],
         'BambusicolaThoracicus_genome.idx': [8900603],
         'CallithrixJacchus_genome.idx': [98299],
         'CentrocercusMinimus_genome.idx': [77819],
         'CentrocercusUrophasianus_genome.idx': [114683],
         'CercopithecusMona_genome.idx': [135163],
         'ChlorocebusSabaeus_genome.idx': [438267],
         'ChrysolophusPictus_genome.idx': [753659],
         'CoturnixJaponica_genome.idx': [126971],
         'CrossoptilonMantchuricum_genome.idx': [163835],
         'CryptoproctaFerox_genome.idx': [26206203],
         'FelisCatus_genome.idx': [319483],
         'GorillaGorilla_genome.idx': [323579],
         'GrammomysDolichurus_genome.idx': [30707707],
         'GrammomysSurdaster_genome.idx': [1282043],
         'GuloGulo_genome.idx': [4288507],
         'HyaenaHyaena_genome.idx': [19730427],
         'HylobatesMoloch_genome.idx': [1023995],
         'HylomyscusAlleni_genome.idx': [35205115],
         'LagopusLeucura_genome.idx': [450555],
         'LagopusMuta_genome.idx': [28667],
         'LophuraNycthemera_genome.idx': [114683],
         'LutraLutra_genome.idx': [20475],
         'LynxRufus_genome.idx': [20475],
         'LyrurusTetrix_genome.idx': [52944891],
         'MacacaMulatta_genome.idx': [126971],
         'MastomysCoucha_genome.idx': [20475],
         'MastomysNatalensis_genome.idx': [19394555],
         'MeleagrisGallopavo_genome.idx': [10027003],
         'MelesMeles_genome.idx': [57339],
         'MiroungaLeonina_genome.idx': [90107],
         'MusCaroli_genome.idx': [192507],
         'MusMinutoides_genome.idx': [4022267],
         'MusMusculus_genome.idx': [0],
         'MusPahari_genome.idx': [163835],
         'MusSpicilegus_genome.idx': [2019323],
         'MusSpretus_genome.idx': [319483],
         'MustelaNigripes_genome.idx': [1269755],
         'NomascusLeucogenys_genome.idx': [147451],
         'OdobenusRosmarus_genome.idx': [212987],
         'PanTroglodytes_genome.idx': [262139],
         'PantheraTigris_genome.idx': [98299],
         'PapioAnubis_genome.idx': [626683],
         'ParadoxurusHermaphroditus_genome.idx': [24698875],
         'PavoCristatus_genome.idx': [65531],
         'PavoMuticus_genome.idx': [167931],
         'PhasianusColchicus_genome.idx': [2134011],
         'PiliocolobusTephrosceles_genome.idx': [2613243],
         'PitheciaPithecia_genome.idx': [50417659],
         'PlecturocebusDonacophilus_genome.idx': [60260347],
         'PongoAbelli_genome.idx': [311291],
         'PraomysDelectorum_genome.idx': [14958587],
         'ProcyonLotor_genome.idx': [3039227],
         'RattusNorvegicus_genome.idx': [69627],
         'RattusRattus_genome.idx': [155643],
         'RhabdomysDilectus_genome.idx': [1880059],
         'RhynchomysSoricoides_genome.idx': [864251],
         'SaimiriBoliviensis_genome.idx': [126971],
         'SpeothosVenaticus_genome.idx': [2867195],
         'SpilogaleGracilis_genome.idx': [23805947],
         'SuricataSuricatta_genome.idx': [3588091],
         'SyrmaticusMikado_genome.idx': [3862523],
         'TrachypithecusFrancoisi_genome.idx': [262139],
         'TympanuchusCupido_genome.idx': [659451],
         'UrsusAmericanus_genome.idx': [6258683],
         'UrsusMaritimus_genome.idx': [249851],
         'VulpesFerrilata_genome.idx': [40955],
         'ZalophusCalifornianus_genome.idx': [20475],
         'DrosophilaSimulans_genome.idx': [35960],
         'DrosophilaMauritiana_genome.idx': [44152],
         'DrosophilaSechellia_genome.idx': [52344],
         'DrosophilaYakuba_genome.idx': [44152],
         'DrosophilaSantomea_genome.idx': [15480],
         'DrosophilaTeissieri_genome.idx': [85112],
         'DrosophilaOrena_genome.idx': [44152],
         'DrosophilaErecta_genome.idx': [52344],
         'DrosophilaEugracilis_genome.idx': [146552],
         'DrosophilaSubpulchrella_genome.idx': [101496],
         'DrosophilaSuzukii_genome.idx': [48248],
         'DrosophilaBiarmipes_genome.idx': [35960],
         'DrosophilaTakahashii_genome.idx': [15480],
         'DrosophilaFicusphila_genome.idx': [68728],
         'DrosophilaCarrolli_genome.idx': [40056],
         'DrosophilaRhopaloa_genome.idx': [31864],
         'DrosophilaKurseongensis_genome.idx': [35960],
         'DrosophilaFuyamai_genome.idx': [72824],
         'DrosophilaElegans_genome.idx': [60536],
         'PanPaniscus_genome.idx': [44152],
         'PongoPygmaeus_genome.idx': [113784],
         'NomascusSiki_genome.idx': [175224],
         'HylobatesPileatus_genome.idx': [64632],
         'SymphalangusSyndactylus_genome.idx': [56440],
         'HoolockLeuconedys_genome.idx': [195704],
         'OtocolobusManul_genome.idx': [0],
         'PrionailurusViverrinus_genome.idx': [0],
         'AcinonyxJubatus_genome.idx': [0],
         'PumaConcolor_genome.idx': [0],
         'NeofelisNebulosa_genome.idx': [0],
         'MungosMungo_genome.idx': [0]}

    expected_index_size = possible_indexes[index_path.split("/")[-1]][0]
    if os.path.getsize(index_path) < expected_index_size:
        return False

    else:
        return True


