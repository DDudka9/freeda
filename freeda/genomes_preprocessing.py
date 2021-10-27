#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 20:57:13 2021

@author: damian
"""


def get_ref_genome_contigs_dict(ref_species):
    """Returns a dictionary of contigs in Genome object (ensembl 104; GRCm39) as keys and
    Sequence-Name or GeneBank-Accn of GCA_000001635.9 mouse genome as values
    
    or
    
    Returns a dictionary of contigs in Genome object (ensembl 104; GRCh38.p13) as keys and
    Sequence-Name or GeneBank-Accn of GCA_000001405.28 human genome as values"""
    
    mouse_names = {"Mm", "mouse", "Mus musculus", "mus musculus"}
    human_names = {"Hs", "human", "Homo sapiens", "homo sapiens"}


    musculus_dict = {'1': 'CM000994.3',
                    '10': 'CM001003.3',
                    '11': 'CM001004.3',
                    '12': 'CM001005.3',
                    '13': 'CM001006.3',
                    '14': 'CM001007.3',
                    '15': 'CM001008.3',
                    '16': 'CM001009.3',
                    '17': 'CM001010.3',
                    '18': 'CM001011.3',
                    '19': 'CM001012.3',
                    '2': 'CM000995.3',
                    '3': 'CM000996.3',
                    '4': 'CM000997.3',
                    '5': 'CM000998.3',
                    '6': 'CM000999.3',
                    '7': 'CM001000.3',
                    '8': 'CM001001.3',
                    '9': 'CM001002.3',
                    'GL456210.1': 'GL456210.1',
                    'GL456211.1': 'GL456211.1',
                    'GL456212.1': 'GL456212.1',
                    'GL456219.1': 'GL456219.1',
                    'GL456221.1': 'GL456221.1',
                    'GL456239.1': 'GL456239.1',
                    'GL456354.1': 'GL456354.1',
                    'GL456372.1': 'GL456372.1',
                    'GL456381.1': 'GL456381.1',
                    'GL456385.1': 'GL456385.1',
                    'JH584295.1': 'JH584295.1',
                    'JH584296.1': 'JH584296.1',
                    'JH584297.1': 'JH584297.1',
                    'JH584298.1': 'JH584298.1',
                    'JH584299.1': 'JH584299.1',
                    'JH584303.1': 'JH584303.1',
                    'JH584304.1': 'JH584304.1',
                    'MT': 'AY172335.1',
                    'X': 'CM001013.3',
                    'Y': 'CM001014.3'}


    sapiens_dict = {'1': 'CM000663.2',
                    '10': 'CM000672.2',
                    '11': 'CM000673.2',
                    '12': 'CM000674.2',
                    '13': 'CM000675.2',
                    '14': 'CM000676.2',
                    '15': 'CM000677.2',
                    '16': 'CM000678.2',
                    '17': 'CM000679.2',
                    '18': 'CM000680.2',
                    '19': 'CM000681.2',
                    '2': 'CM000664.2',
                    '20': 'CM000682.2',
                    '21': 'CM000683.2',
                    '22': 'CM000684.2',
                    '3': 'CM000665.2',
                    '4': 'CM000666.2',
                    '5': 'CM000667.2',
                    '6': 'CM000668.2',
                    '7': 'CM000669.2',
                    '8': 'CM000670.2',
                    '9': 'CM000671.2',
                    'GL000009.2': 'GL000009.2',
                    'GL000194.1': 'GL000194.1',
                    'GL000195.1': 'GL000195.1',
                    'GL000205.2': 'GL000205.2',
                    'GL000213.1': 'GL000213.1',
                    'GL000216.2': 'GL000216.2',
                    'GL000218.1': 'GL000218.1',
                    'GL000219.1': 'GL000219.1',
                    'GL000220.1': 'GL000220.1',
                    'GL000225.1': 'GL000225.1',
                    'KI270442.1': 'KI270442.1',
                    'KI270711.1': 'KI270711.1',
                    'KI270713.1': 'KI270713.1',
                    'KI270721.1': 'KI270721.1',
                    'KI270726.1': 'KI270726.1',
                    'KI270727.1': 'KI270727.1',
                    'KI270728.1': 'KI270728.1',
                    'KI270731.1': 'KI270731.1',
                    'KI270733.1': 'KI270733.1',
                    'KI270734.1': 'KI270734.1',
                    'KI270744.1': 'KI270744.1',
                    'KI270750.1': 'KI270750.1',
                    'MT': 'J01415.2',
                    'X': 'CM000685.2',
                    'Y': 'CM000686.2'}

    # default is mouse
    if ref_species in mouse_names:
        ref_genome_contigs_dict = musculus_dict
    
    elif ref_species in human_names:
        ref_genome_contigs_dict = sapiens_dict
        
    return ref_genome_contigs_dict


def get_names(ref_species, ref_genome=False):
    """Gets species, genomes names and accession numbers used for FREEDA analysis"""

    mouse_dict = {"Mm": (#("Mi", "SPICILEGUS_genome", "GCA_003336285.1"),
                      #("Ms", "SPRETUS_genome", "GCA_001624865.1"),
                      #("Mc", "CAROLI_genome", "GCA_900094665.2"),
                      #("Mu", "MINUTOIDES_genome", "GCA_902729485.2"),
                      #("Mp", "PAHARI_genome", "GCA_900095145.2"),
                      #("Ay", "SYLVATICUS_genome", "GCA_001305905.1"),
                      #("Ap", "SPECIOSUS_genome", "GCA_002335545.1"),
                      # ("Ha", "ALLENI_genome", "GCA_019843855.1"),  # got rid of this genome -> too fragmented
                      #("Pd", "DELECTORUM_genome", "GCA_019843815.1"),
                      #("Mn", "NATALENSIS_genome", "GCA_019843795.1"),
                      #("Mo", "COUCHA_genome", "GCA_008632895.1"),
                      #("Gd", "DOLICHURUS_genome", "GCA_019843835.1"),  # this one is almost as bad as Ha
                      #("Gs", "SURDASTER_genome", "GCA_004785775.1"),
                      #("An", "ARVICANTHIS_genome", "GCA_011762505.1"),
                      #("Rd", "DILECTUS_genome", "GCA_019844195.1"),
                      #("Rs", "SORICOIDES_genome", "GCA_019843965.1"),
                      ("Rr", "RATTUS_genome", "GCA_011800105.1"),  # Rattus rattus (Black rat)
                      ("Rn", "NORVEGICUS_genome", "GCA_015227675.2"))}  # Rattus norvegicus (Norway rat)

    human_dict = {"Hs": (("Pt", "TROGLODYTES_genome", "GCA_002880755.3"),
                      ("Gg", "GORILLA_genome", "GCA_008122165.1"),
                      ("Pb", "ABELII_genome", "GCA_002880775.3"),
                      ("Ne", "LEUCOGENYS_genome", "GCA_006542625.1"),
                      ("Mu", "MULATTA_genome", "GCA_008058575.1"),
                      ("Pu", "ANUBIS_genome", "GCA_008728515.1"),
                      ("Cs", "SABAEUS_genome", "GCA_015252025.1"),
                      ("Cj", "JACCHUS_genome", "GCA_011100535.2"))}

    mouse_ref_dict = {"Mm": (("Mm", "MUSCULUS_genome", "GCA_000001635.9"))} # GeneBank GRCm39
    human_ref_dict = {"Hs": (("Hs", "SAPIENS_genome", "GCA_000001405.28"))} # RefSeq GCF_000001405.39

    if ref_species == "Mm" and ref_genome is False:
        genomes_dict = mouse_dict
    elif ref_species == "Mm" and ref_genome is True:
        genomes_dict = mouse_ref_dict
    elif ref_species == "Hs" and ref_genome is False:
        genomes_dict = human_dict
    elif ref_species == "Hs" and ref_genome is True:
        genomes_dict = human_ref_dict
    else:
        print("Something went wrong")

    # collect genomes
    all_genomes = []
    for ref, species in genomes_dict.items():
        for genome in species: # genome = ("Mi", "SPICILEGUS_genome", "GCA_003336285.1")
            all_genomes.append(genome)

    return all_genomes










"""

def map_assembly_contigs(wdir):

    import pyensembl

    release = 104
    species = "mus musculus"

    ensembl = pyensembl.EnsemblRelease(release, species)

    # get contigs available in ensembl release
    contigs = ensembl.contigs()

    # get all available contigs in the NCBI assembly
    with open(wdir + "MUSCULUS_genome_contigs.txt", "r") as f:
        rows = (line.split('\t') for line in f)
        # need index 3 from each list
        all_contigs = {row[0]: row[1:] for row in rows}

    # get contigs from indexed reference assembly
    with open(wdir + "Reference_genomes/MUSCULUS_genome.fasta.fai", "r") as f:
        rows = (line.split('\t') for line in f)
        # need keys
        indexed_contigs = {row[0]: row[1:] for row in rows}

    mapped_contigs_dict = {}
    for contig in contigs:
        for indexed_contig in indexed_contigs:
            for available_contig, features in all_contigs.items():
                if available_contig == contig and indexed_contig == features[3]:
                    mapped_contigs_dict[contig] = indexed_contig
                elif contig == features[3]:
                    mapped_contigs_dict[contig] = features[3]

    return mapped_contigs_dict

# deprecated (GRCm38.6)
    
    musculus_dict_RefSeq = {'1': 'NC_000067.6',
                '10': 'NC_000076.6',
                '11': 'NC_000077.6',
                '12': 'NC_000078.6',
                '13': 'NC_000079.6',
                '14': 'NC_000080.6',
                '15': 'NC_000081.6',
                '16': 'NC_000082.6',
                '17': 'NC_000083.6',
                '18': 'NC_000084.6',
                '19': 'NC_000085.6',
                '2': 'NC_000068.7',
                '3': 'NC_000069.6',
                '4': 'NC_000070.6',
                '5': 'NC_000071.6',
                '6': 'NC_000072.6',
                '7': 'NC_000073.6',
                '8': 'NC_000074.6',
                '9': 'NC_000075.6',
                'GL456210.1': 'NT_166280.1',
                'GL456211.1': 'NT_166281.1',
                'GL456212.1': 'NT_166282.1',
                'GL456216.1': 'NT_166291.1',
                'GL456219.1': 'NT_166307.1',
                'GL456221.1': 'NT_162750.1',
                'GL456233.1': 'NT_165789.2',
                'GL456239.1': 'NT_166338.1',
                'GL456350.1': 'NT_166434.1',
                'GL456354.1': 'NT_166438.1',
                'GL456372.1': 'NT_166456.1',
                'GL456381.1': 'NT_166465.1',
                'GL456385.1': 'NT_166469.1',
                'JH584292.1': 'NT_187052.1',
                'JH584293.1': 'NT_187053.1',
                'JH584294.1': 'NT_187054.1',
                'JH584295.1': 'NT_187055.1',
                'JH584296.1': 'NT_187056.1',
                'JH584297.1': 'NT_187057.1',
                'JH584298.1': 'NT_187058.1',
                'JH584299.1': 'NT_187059.1',
                'JH584303.1': 'NT_187063.1',
                'JH584304.1': 'NT_187064.1',
                'MT': 'NC_005089.1',
                'X': 'NC_000086.7',
                'Y': 'NC_000087.7'}



    musculus_dict_GenBank = {'1': 'CM000994.2',
                     '10': 'CM001003.2',
                     '11': 'CM001004.2',
                     '12': 'CM001005.2',
                     '13': 'CM001006.2',
                     '14': 'CM001007.2',
                     '15': 'CM001008.2',
                     '16': 'CM001009.2',
                     '17': 'CM001010.2',
                     '18': 'CM001011.2',
                     '19': 'CM001012.2',
                     '2': 'CM000995.2',
                     '3': 'CM000996.2',
                     '4': 'CM000997.2',
                     '5': 'CM000998.2',
                     '6': 'CM000999.2',
                     '7': 'CM001000.2',
                     '8': 'CM001001.2',
                     '9': 'CM001002.2',
                     'GL456210.1': 'GL456210.1',
                     'GL456211.1': 'GL456211.1',
                     'GL456212.1': 'GL456212.1',
                     'GL456216.1': 'GL456216.1',
                     'GL456219.1': 'GL456219.1',
                     'GL456221.1': 'GL456221.1',
                     'GL456233.1': 'GL456233.1',
                     'GL456239.1': 'GL456239.1',
                     'GL456350.1': 'GL456350.1',
                     'GL456354.1': 'GL456354.1',
                     'GL456372.1': 'GL456372.1',
                     'GL456381.1': 'GL456381.1',
                     'GL456385.1': 'GL456385.1',
                     'JH584292.1': 'JH584292.1',
                     'JH584293.1': 'JH584293.1',
                     'JH584294.1': 'JH584294.1',
                     'JH584295.1': 'JH584295.1',
                     'JH584296.1': 'JH584296.1',
                     'JH584297.1': 'JH584297.1',
                     'JH584298.1': 'JH584298.1',
                     'JH584299.1': 'JH584299.1',
                     'JH584303.1': 'JH584303.1',
                     'JH584304.1': 'JH584304.1',
                     'MT': 'AY172335.1',
                     'X': 'CM001013.2',
                     'Y': 'CM001014.2'}


Notes from meeting with NCBI Datasets developpers on 09/24/2021


    unzip - p
    dataset.zip ‘chr *.fna
    ' > all_chr_files.fna

    Correction(hopefully): unzip - p
    ncbi_dataset.zip
    '*/chr*.fna' > all_chr_files.fna

    https: // anaconda.org / conda - forge / ncbi - datasets - cli

    datasets
    summary
    gene
    symbol
    sumo1 | jq. | less

"""
