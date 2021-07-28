#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 20:57:13 2021

@author: damian
"""

musculus_dict = {'1': 'CM000994.2',
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

def get_reference_genome_contigs_dict(original_species):
    """Returns a dictionary of contigs in Genome object (ensembl 100; GRCm38.p6) as keys and 
    Sequence-Name or GeneBank-Accn of GCA_000001635.8 mouse genome as values
    
    or
    
    Returns a dictionary of contigs in Genome object (ensembl 100; GRCh38.p13) as keys and 
    Sequence-Name or GeneBank-Accn of GCA_000001405.28 human genome as values"""
    
    mouse_names = {"Mm", "mouse", "Mus musculus", "mus musculus"}
    human_names = {"Hs", "human", "Homo sapiens", "homo sapiens"}
    
    # default is mouse
    if original_species in mouse_names:
        reference_genome_contigs_dict = musculus_dict
    
    elif original_species in human_names:
        reference_genome_contigs_dict = sapiens_dict
        
    return reference_genome_contigs_dict
























