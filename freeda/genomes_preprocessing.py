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


    mouse_dict = {'1': 'CM000994.3',
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

    rat_dict = {'1': 'CM000072.5',
                 '10': 'CM000081.5',
                 '11': 'CM000082.5',
                 '12': 'CM000083.5',
                 '13': 'CM000084.5',
                 '14': 'CM000085.5',
                 '15': 'CM000086.5',
                 '16': 'CM000087.5',
                 '17': 'CM000088.5',
                 '18': 'CM000089.5',
                 '19': 'CM000090.5',
                 '2': 'CM000073.5',
                 '20': 'CM000091.5',
                 '3': 'CM000074.5',
                 '4': 'CM000075.5',
                 '5': 'CM000076.5',
                 '6': 'CM000077.5',
                 '7': 'CM000078.5',
                 '8': 'CM000079.5',
                 '9': 'CM000080.5',
                 'AABR07022258.1': 'AABR07022258.1',
                 'AABR07022620.1': 'AABR07022620.1',
                 'AABR07022926.1': 'AABR07022926.1',
                 'AABR07024031.1': 'AABR07024031.1',
                 'AABR07024032.1': 'AABR07024032.1',
                 'AABR07024040.1': 'AABR07024040.1',
                 'AABR07024041.1': 'AABR07024041.1',
                 'AABR07024044.1': 'AABR07024044.1',
                 'AABR07024102.1': 'AABR07024102.1',
                 'AABR07024104.1': 'AABR07024104.1',
                 'AABR07024106.1': 'AABR07024106.1',
                 'AABR07024115.1': 'AABR07024115.1',
                 'AABR07024119.1': 'AABR07024119.1',
                 'AABR07024120.1': 'AABR07024120.1',
                 'AABR07024122.1': 'AABR07024122.1',
                 'AABR07024124.1': 'AABR07024124.1',
                 'AABR07024203.1': 'AABR07024203.1',
                 'AABR07024206.1': 'AABR07024206.1',
                 'AABR07024263.1': 'AABR07024263.1',
                 'AABR07024264.1': 'AABR07024264.1',
                 'AABR07024291.1': 'AABR07024291.1',
                 'AABR07024382.1': 'AABR07024382.1',
                 'AABR07024399.1': 'AABR07024399.1',
                 'AABR07024421.1': 'AABR07024421.1',
                 'AABR07046136.1': 'AABR07046136.1',
                 'AABR07046231.1': 'AABR07046231.1',
                 'AABR07046248.1': 'AABR07046248.1',
                 'AABR07046563.1': 'AABR07046563.1',
                 'KL567891.1': 'KL567891.1',
                 'KL567892.1': 'KL567892.1',
                 'KL567906.1': 'KL567906.1',
                 'KL567908.1': 'KL567908.1',
                 'KL567939.1': 'KL567939.1',
                 'KL567971.1': 'KL567971.1',
                 'KL567988.1': 'KL567988.1',
                 'KL568001.1': 'KL568001.1',
                 'KL568002.1': 'KL568002.1',
                 'KL568009.1': 'KL568009.1',
                 'KL568013.1': 'KL568013.1',
                 'KL568017.1': 'KL568017.1',
                 'KL568024.1': 'KL568024.1',
                 'KL568037.1': 'KL568037.1',
                 'KL568076.1': 'KL568076.1',
                 'KL568084.1': 'KL568084.1',
                 'KL568092.1': 'KL568092.1',
                 'KL568097.1': 'KL568097.1',
                 'KL568119.1': 'KL568119.1',
                 'KL568122.1': 'KL568122.1',
                 'KL568125.1': 'KL568125.1',
                 'KL568128.1': 'KL568128.1',
                 'KL568132.1': 'KL568132.1',
                 'KL568139.1': 'KL568139.1',
                 'KL568141.1': 'KL568141.1',
                 'KL568142.1': 'KL568142.1',
                 'KL568143.1': 'KL568143.1',
                 'KL568146.1': 'KL568146.1',
                 'KL568147.1': 'KL568147.1',
                 'KL568148.1': 'KL568148.1',
                 'KL568149.1': 'KL568149.1',
                 'KL568150.1': 'KL568150.1',
                 'KL568151.1': 'KL568151.1',
                 'KL568152.1': 'KL568152.1',
                 'KL568156.1': 'KL568156.1',
                 'KL568157.1': 'KL568157.1',
                 'KL568158.1': 'KL568158.1',
                 'KL568159.1': 'KL568159.1',
                 'KL568160.1': 'KL568160.1',
                 'KL568161.1': 'KL568161.1',
                 'KL568162.1': 'KL568162.1',
                 'KL568169.1': 'KL568169.1',
                 'KL568174.1': 'KL568174.1',
                 'KL568194.1': 'KL568194.1',
                 'KL568195.1': 'KL568195.1',
                 'KL568198.1': 'KL568198.1',
                 'KL568199.1': 'KL568199.1',
                 'KL568205.1': 'KL568205.1',
                 'KL568208.1': 'KL568208.1',
                 'KL568212.1': 'KL568212.1',
                 'KL568217.1': 'KL568217.1',
                 'KL568218.1': 'KL568218.1',
                 'KL568221.1': 'KL568221.1',
                 'KL568231.1': 'KL568231.1',
                 'KL568233.1': 'KL568233.1',
                 'KL568234.1': 'KL568234.1',
                 'KL568242.1': 'KL568242.1',
                 'KL568244.1': 'KL568244.1',
                 'KL568257.1': 'KL568257.1',
                 'KL568263.1': 'KL568263.1',
                 'KL568274.1': 'KL568274.1',
                 'KL568281.1': 'KL568281.1',
                 'KL568295.1': 'KL568295.1',
                 'KL568297.1': 'KL568297.1',
                 'KL568300.1': 'KL568300.1',
                 'KL568304.1': 'KL568304.1',
                 'KL568305.1': 'KL568305.1',
                 'KL568306.1': 'KL568306.1',
                 'KL568307.1': 'KL568307.1',
                 'KL568325.1': 'KL568325.1',
                 'KL568337.1': 'KL568337.1',
                 'KL568367.1': 'KL568367.1',
                 'KL568371.1': 'KL568371.1',
                 'KL568374.1': 'KL568374.1',
                 'KL568381.1': 'KL568381.1',
                 'KL568401.1': 'KL568401.1',
                 'KL568405.1': 'KL568405.1',
                 'KL568409.1': 'KL568409.1',
                 'KL568410.1': 'KL568410.1',
                 'KL568411.1': 'KL568411.1',
                 'KL568414.1': 'KL568414.1',
                 'KL568417.1': 'KL568417.1',
                 'KL568418.1': 'KL568418.1',
                 'KL568423.1': 'KL568423.1',
                 'KL568430.1': 'KL568430.1',
                 'KL568431.1': 'KL568431.1',
                 'KL568432.1': 'KL568432.1',
                 'KL568435.1': 'KL568435.1',
                 'KL568438.1': 'KL568438.1',
                 'KL568439.1': 'KL568439.1',
                 'KL568447.1': 'KL568447.1',
                 'KL568449.1': 'KL568449.1',
                 'KL568451.1': 'KL568451.1',
                 'KL568458.1': 'KL568458.1',
                 'KL568460.1': 'KL568460.1',
                 'KL568463.1': 'KL568463.1',
                 'KL568466.1': 'KL568466.1',
                 'KL568468.1': 'KL568468.1',
                 'KL568470.1': 'KL568470.1',
                 'KL568472.1': 'KL568472.1',
                 'KL568473.1': 'KL568473.1',
                 'KL568483.1': 'KL568483.1',
                 'KL568487.1': 'KL568487.1',
                 'KL568488.1': 'KL568488.1',
                 'KL568489.1': 'KL568489.1',
                 'KL568496.1': 'KL568496.1',
                 'KL568497.1': 'KL568497.1',
                 'KL568507.1': 'KL568507.1',
                 'KL568511.1': 'KL568511.1',
                 'KL568512.1': 'KL568512.1',
                 'KL568518.1': 'KL568518.1',
                 'MT': 'AY172581.1',
                 'X': 'CM000092.5',
                 'Y': 'CM002824.1'}

    human_dict = {'1': 'CM000663.2',
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

    cat_dict = {'A1': 'NC_018723.1',
     'A2': 'NC_018724.1',
     'A3': 'NC_018725.1',
     'B1': 'NC_018726.1',
     'B2': 'NC_018727.1',
     'B3': 'NC_018728.1',
     'B4': 'NC_018729.1',
     'C1': 'NC_018730.1',
     'C2': 'NC_018731.1',
     'D1': 'NC_018732.1',
     'D2': 'NC_018733.1',
     'D3': 'NC_018734.1',
     'D4': 'NC_018735.1',
     'E1': 'NC_018736.1',
     'E2': 'NC_018737.1',
     'E3': 'NC_018738.1',
     'F1': 'NC_018739.1',
     'F2': 'NC_018740.1',
     'JH408314.1': 'NW_004065384.1',
     'JH408315.1': 'NW_004065385.1',
     'JH408316.1': 'NW_004065386.1',
     'JH408317.1': 'NW_004065387.1',
     'JH408319.1': 'NW_004065389.1',
     'JH408320.1': 'NW_004065390.1',
     'JH408321.1': 'NW_004065391.1',
     'JH408329.1': 'NW_004065399.1',
     'JH408367.1': 'NW_004065437.1',
     'JH408372.1': 'NW_004065442.1',
     'JH408373.1': 'NW_004065443.1',
     'JH408375.1': 'NW_004065445.1',
     'JH408376.1': 'NW_004065446.1',
     'JH408378.1': 'NW_004065448.1',
     'JH408380.1': 'NW_004065450.1',
     'JH408408.1': 'NW_004065478.1',
     'JH408414.1': 'NW_004065484.1',
     'JH408425.1': 'NW_004065495.1',
     'JH408484.1': 'NW_004065554.1',
     'JH408486.1': 'NW_004065556.1',
     'JH408493.1': 'NW_004065563.1',
     'JH408514.1': 'NW_004065584.1',
     'JH408522.1': 'NW_004065592.1',
     'JH408533.1': 'NW_004065603.1',
     'JH408547.1': 'NW_004065617.1',
     'JH408548.1': 'NW_004065618.1',
     'JH408551.1': 'NW_004065621.1',
     'JH408553.1': 'NW_004065623.1',
     'JH408554.1': 'NW_004065624.1',
     'JH408571.1': 'NW_004065641.1',
     'JH408572.1': 'NW_004065642.1',
     'JH408577.1': 'NW_004065647.1',
     'JH408579.1': 'NW_004065649.1',
     'JH408580.1': 'NW_004065650.1',
     'JH408582.1': 'NW_004065652.1',
     'JH408593.1': 'NW_004065663.1',
     'JH408601.1': 'NW_004065671.1',
     'JH408634.1': 'NW_004065704.1',
     'JH408635.1': 'NW_004065705.1',
     'JH408636.1': 'NW_004065706.1',
     'JH408642.1': 'NW_004065712.1',
     'JH408646.1': 'NW_004065716.1',
     'JH408653.1': 'NW_004065723.1',
     'JH408672.1': 'NW_004065742.1',
     'JH408674.1': 'NW_004065744.1',
     'JH408682.1': 'NW_004065752.1',
     'JH408685.1': 'NW_004065755.1',
     'JH408686.1': 'NW_004065756.1',
     'JH408688.1': 'NW_004065758.1',
     'JH408690.1': 'NW_004065760.1',
     'JH408714.1': 'NW_004065784.1',
     'JH408723.1': 'NW_004065793.1',
     'JH408730.1': 'NW_004065800.1',
     'JH408731.1': 'NW_004065801.1',
     'JH408732.1': 'NW_004065802.1',
     'JH408733.1': 'NW_004065803.1',
     'JH408737.1': 'NW_004065807.1',
     'JH408758.1': 'NW_004065828.1',
     'JH408768.1': 'NW_004065838.1',
     'JH408779.1': 'NW_004065849.1',
     'JH408785.1': 'NW_004065855.1',
     'JH408787.1': 'NW_004065857.1',
     'JH408788.1': 'NW_004065858.1',
     'JH408801.1': 'NW_004065871.1',
     'JH408808.1': 'NW_004065878.1',
     'JH408818.1': 'NW_004065888.1',
     'JH408842.1': 'NW_004065912.1',
     'JH408843.1': 'NW_004065913.1',
     'JH408869.1': 'NW_004065939.1',
     'JH408870.1': 'NW_004065940.1',
     'JH408873.1': 'NW_004065943.1',
     'JH408877.1': 'NW_004065947.1',
     'JH408889.1': 'NW_004065959.1',
     'JH408897.1': 'NW_004065967.1',
     'JH408900.1': 'NW_004065970.1',
     'JH408902.1': 'NW_004065972.1',
     'JH408903.1': 'NW_004065973.1',
     'JH408904.1': 'NW_004065974.1',
     'JH408905.1': 'NW_004065975.1',
     'JH408907.1': 'NW_004065977.1',
     'JH408929.1': 'NW_004065999.1',
     'JH408930.1': 'NW_004066000.1',
     'JH408931.1': 'NW_004066001.1',
     'JH408935.1': 'NW_004066005.1',
     'JH408941.1': 'NW_004066011.1',
     'JH408942.1': 'NW_004066012.1',
     'JH408946.1': 'NW_004066016.1',
     'JH408948.1': 'NW_004066018.1',
     'JH408957.1': 'NW_004066027.1',
     'JH408959.1': 'NW_004066029.1',
     'JH408960.1': 'NW_004066030.1',
     'JH408961.1': 'NW_004066031.1',
     'JH408963.1': 'NW_004066033.1',
     'JH408972.1': 'NW_004066042.1',
     'JH408973.1': 'NW_004066043.1',
     'JH408974.1': 'NW_004066044.1',
     'JH408975.1': 'NW_004066045.1',
     'JH408976.1': 'NW_004066046.1',
     'JH408977.1': 'NW_004066047.1',
     'JH408987.1': 'NW_004066057.1',
     'JH408988.1': 'NW_004066058.1',
     'JH408995.1': 'NW_004066065.1',
     'JH409009.1': 'NW_004066079.1',
     'JH409016.1': 'NW_004066086.1',
     'JH409018.1': 'NW_004066088.1',
     'JH409020.1': 'NW_004066090.1',
     'JH409044.1': 'NW_004066114.1',
     'JH409055.1': 'NW_004066125.1',
     'JH409096.1': 'NW_004066166.1',
     'JH409109.1': 'NW_004066179.1',
     'JH409123.1': 'NW_004066193.1',
     'JH409232.1': 'NW_004066302.1',
     'JH409246.1': 'NW_004066316.1',
     'JH409275.1': 'NW_004066345.1',
     'JH409321.1': 'NW_004066391.1',
     'JH409365.1': 'NW_004066435.1',
     'JH409408.1': 'NW_004066478.1',
     'JH409432.1': 'NW_004066502.1',
     'JH409652.1': 'NW_004066722.1',
     'JH409676.1': 'NW_004066746.1',
     'JH409756.1': 'NW_004066826.1',
     'JH409757.1': 'NW_004066827.1',
     'JH409775.1': 'NW_004066845.1',
     'JH409781.1': 'NW_004066851.1',
     'JH409793.1': 'NW_004066863.1',
     'JH409826.1': 'NW_004066896.1',
     'JH409937.1': 'NW_004067007.1',
     'JH409952.1': 'NW_004067022.1',
     'JH410053.1': 'NW_004067123.1',
     'JH410060.1': 'NW_004067130.1',
     'JH410159.1': 'NW_004067229.1',
     'JH410221.1': 'NW_004067291.1',
     'JH410393.1': 'NW_004067463.1',
     'JH410469.1': 'NW_004067539.1',
     'JH410575.1': 'NW_004067645.1',
     'JH410582.1': 'NW_004067652.1',
     'JH410585.1': 'NW_004067655.1',
     'JH410651.1': 'NW_004067721.1',
     'JH410698.1': 'NW_004067768.1',
     'JH410715.1': 'NW_004067785.1',
     'JH410722.1': 'NW_004067792.1',
     'JH410725.1': 'NW_004067795.1',
     'JH410871.1': 'NW_004067941.1',
     'JH410880.1': 'NW_004067950.1',
     'JH410902.1': 'NW_004067972.1',
     'JH410918.1': 'NW_004067988.1',
     'JH410921.1': 'NW_004067991.1',
     'JH410946.1': 'NW_004068016.1',
     'JH410965.1': 'NW_004068035.1',
     'JH410968.1': 'NW_004068038.1',
     'JH411006.1': 'NW_004068076.1',
     'JH411061.1': 'NW_004068131.1',
     'JH411092.1': 'NW_004068162.1',
     'JH411094.1': 'NW_004068164.1',
     'JH411117.1': 'NW_004068187.1',
     'JH411125.1': 'NW_004068195.1',
     'JH411136.1': 'NW_004068206.1',
     'JH411139.1': 'NW_004068209.1',
     'JH411178.1': 'NW_004068248.1',
     'JH411188.1': 'NW_004068258.1',
     'JH411233.1': 'NW_004068303.1',
     'JH411285.1': 'NW_004068355.1',
     'JH411297.1': 'NW_004068367.1',
     'JH411298.1': 'NW_004068368.1',
     'JH411318.1': 'NW_004068388.1',
     'JH411322.1': 'NW_004068392.1',
     'JH411389.1': 'NW_004068459.1',
     'JH411489.1': 'NW_004068559.1',
     'JH411490.1': 'NW_004068560.1',
     'JH411491.1': 'NW_004068561.1',
     'JH411525.1': 'NW_004068595.1',
     'JH411592.1': 'NW_004068662.1',
     'JH411605.1': 'NW_004068675.1',
     'JH411664.1': 'NW_004068734.1',
     'JH411748.1': 'NW_004068818.1',
     'JH411754.1': 'NW_004068824.1',
     'JH411779.1': 'NW_004068849.1',
     'JH411790.1': 'NW_004068860.1',
     'JH411805.1': 'NW_004068875.1',
     'JH411853.1': 'NW_004068923.1',
     'JH411911.1': 'NW_004068981.1',
     'JH411917.1': 'NW_004068987.1',
     'JH411936.1': 'NW_004069006.1',
     'JH411971.1': 'NW_004069041.1',
     'JH412049.1': 'NW_004069119.1',
     'JH412128.1': 'NW_004069198.1',
     'JH412145.1': 'NW_004069215.1',
     'JH412174.1': 'NW_004069244.1',
     'JH412201.1': 'NW_004069271.1',
     'JH412249.1': 'NW_004069319.1',
     'JH412414.1': 'NW_004069484.1',
     'JH412471.1': 'NW_004069541.1',
     'JH412472.1': 'NW_004069542.1',
     'JH412473.1': 'NW_004069543.1',
     'JH412484.1': 'NW_004069554.1',
     'JH412486.1': 'NW_004069556.1',
     'JH412491.1': 'NW_004069561.1',
     'JH412495.1': 'NW_004069565.1',
     'JH412514.1': 'NW_004069584.1',
     'JH412522.1': 'NW_004069592.1',
     'JH412547.1': 'NW_004069617.1',
     'JH412564.1': 'NW_004069634.1',
     'JH412587.1': 'NW_004069657.1',
     'JH412600.1': 'NW_004069670.1',
     'JH412602.1': 'NW_004069672.1',
     'JH412607.1': 'NW_004069677.1',
     'JH412614.1': 'NW_004069684.1',
     'JH412624.1': 'NW_004069694.1',
     'JH412625.1': 'NW_004069695.1',
     'JH412631.1': 'NW_004069701.1',
     'JH412632.1': 'NW_004069702.1',
     'JH412637.1': 'NW_004069707.1',
     'JH412642.1': 'NW_004069712.1',
     'JH412643.1': 'NW_004069713.1',
     'JH412646.1': 'NW_004069716.1',
     'JH412654.1': 'NW_004069724.1',
     'JH412659.1': 'NW_004069729.1',
     'JH412703.1': 'NW_004069773.1',
     'JH412707.1': 'NW_004069777.1',
     'JH412728.1': 'NW_004069798.1',
     'JH412731.1': 'NW_004069801.1',
     'JH412748.1': 'NW_004069818.1',
     'JH412782.1': 'NW_004069852.1',
     'JH412808.1': 'NW_004069878.1',
     'JH412987.1': 'NW_004070057.1',
     'JH413057.1': 'NW_004070127.1',
     'JH413374.1': 'NW_004070444.1',
     'JH413381.1': 'NW_004070451.1',
     'JH413457.1': 'NW_004070527.1',
     'JH413521.1': 'NW_004070591.1',
     'JH413580.1': 'NW_004070650.1',
     'JH413581.1': 'NW_004070651.1',
     'JH413642.1': 'NW_004070712.1',
     'JH413673.1': 'NW_004070743.1',
     'JH413677.1': 'NW_004070747.1',
     'JH413683.1': 'NW_004070753.1',
     'JH413737.1': 'NW_004070807.1',
     'JH413764.1': 'NW_004070834.1',
     'MT': 'NC_001700.1',
     'X': 'NC_018741.1'}

    # mock
    dog_dict = {'1': 'CM000994.3',
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


    # default is mouse
    if ref_species == "Mm":
        ref_genome_contigs_dict = mouse_dict

    elif ref_species == "Rn":
        ref_genome_contigs_dict = rat_dict
    
    elif ref_species == "Hs":
        ref_genome_contigs_dict = human_dict

    elif ref_species == "Fc":
        ref_genome_contigs_dict = cat_dict

    elif ref_species == "Cf":
        ref_genome_contigs_dict = dog_dict
        
    return ref_genome_contigs_dict


def get_names(ref_species, ref_genome=False):
    """Gets species, genomes names and accession numbers used for FREEDA analysis"""

    mouse_dict = {"Mm": (("Mi", "SPICILEGUS_genome", "GCA_003336285.1"),
                      ("Ms", "SPRETUS_genome", "GCA_001624865.1"),
                      ("Mc", "CAROLI_genome", "GCA_900094665.2"),
                      ("Mu", "MINUTOIDES_genome", "GCA_902729485.2"),
                      ("Mp", "PAHARI_genome", "GCA_900095145.2"),
                      ("Ay", "SYLVATICUS_genome", "GCA_001305905.1"),
                      ("Ap", "SPECIOSUS_genome", "GCA_002335545.1"),
                      ("Ha", "ALLENI_genome", "GCA_019843855.1"),  # got rid of this genome -> too fragmented
                      ("Pd", "DELECTORUM_genome", "GCA_019843815.1"),
                      ("Mn", "NATALENSIS_genome", "GCA_019843795.1"),
                      ("Mo", "COUCHA_genome", "GCA_008632895.1"),
                      ("Gd", "DOLICHURUS_genome", "GCA_019843835.1"),  # this one is almost as bad as Ha
                      ("Gs", "SURDASTER_genome", "GCA_004785775.1"),
                      ("An", "ARVICANTHIS_genome", "GCA_011762505.1"),
                      ("Rd", "DILECTUS_genome", "GCA_019844195.1"),
                      ("Rs", "SORICOIDES_genome", "GCA_019843965.1"),
                      ("Rr", "RATTUS_genome", "GCA_011064425.1"),  # Rattus rattus (Black rat)  alt -> GCA_011800105.1
                      ("Rn", "NORVEGICUS_genome", "GCA_000001895.4"))}  # Previously used : GCA_015227675.2

    rat_dict = {"Rn": (("Rr", "RATTUS_genome", "GCA_011064425.1"),  # Rattus rattus (Black rat)  alt -> GCA_011800105.1
                         ("Rs", "SORICOIDES_genome", "GCA_019843965.1"),
                         ("Rd", "DILECTUS_genome", "GCA_019844195.1"),
                         ("An", "ARVICANTHIS_genome", "GCA_011762505.1"),
                         ("Gs", "SURDASTER_genome", "GCA_004785775.1"),
                         ("Gd", "DOLICHURUS_genome", "GCA_019843835.1"),  # this one is almost as bad as Ha
                         ("Mo", "COUCHA_genome", "GCA_008632895.1"),
                         ("Mn", "NATALENSIS_genome", "GCA_019843795.1"),
                         ("Pd", "DELECTORUM_genome", "GCA_019843815.1"),
                         ("Ha", "ALLENI_genome", "GCA_019843855.1"),  # got rid of this genome -> too fragmented
                         ("Ap", "SPECIOSUS_genome", "GCA_002335545.1"),
                         ("Ay", "SYLVATICUS_genome", "GCA_001305905.1"),
                         ("Mp", "PAHARI_genome", "GCA_900095145.2"),
                         ("Mu", "MINUTOIDES_genome", "GCA_902729485.2"),
                         ("Mc", "CAROLI_genome", "GCA_900094665.2"),
                         ("Ms", "SPRETUS_genome", "GCA_001624865.1"),
                         ("Mi", "SPICILEGUS_genome", "GCA_003336285.1"),
                         ("Mm", "MUSCULUS_genome", "GCA_000001635.9"))}

    human_dict = {"Hs": (("Pt", "TROGLODYTES_genome", "GCA_002880755.3"),
                      ("Gg", "GORILLA_genome", "GCA_008122165.1"),
                      ("Pb", "ABELII_genome", "GCA_002880775.3"),
                      ("Ne", "LEUCOGENYS_genome", "GCA_006542625.1"),
                      ("Mu", "MULATTA_genome", "GCA_008058575.1"),
                      ("Pu", "ANUBIS_genome", "GCA_008728515.1"),
                      ("Pp", "PITHECIA_genome", "GCA_004026645.1"),   # added White-faced saki
                      ("Cs", "SABAEUS_genome", "GCA_015252025.1"),
                      ("Pi", "TEPHROSCELES_genome", "GCA_002776525.3"),  # added Colobus monkey
                      ("An", "NANCYMAAE_genome", "GCA_000952055.2"),  # added Ma's night monkey
                      ("Cj", "JACCHUS_genome", "GCA_011100535.2"),
                      ("Sb", "BOLIVIENSIS_genome", "GCA_016699345.1"),  # added Squirrel monkey
                      ("Ag", "GEOFFROYI_genome", "GCA_004024785.1"))}  # added Spider monkey
                      # ("Mm", "MURINUS_genome", "GCA_000165445.3"),  # added Mouse lemur
                      # ("Og", "GARNETTI_genome", "GCA_000181295.3"))}   # added Galago lemur

    # Felidae (diverged 15 myo - too narrow)
    cat_dict = {"Fc": (("Fh", "Felischaus_genome", "GCA_019924945.1"),  # (jungle cat)
                      ("Pt", "Pantheratigris_genome", "GCA_018350195.2"),  # (tiger)
                      ("Pl", "Pantheraleo_genome", "GCA_018350215.1"),   # (lion)
                      ("Pp", "Pantherapardus_genome", "GCA_001857705.1"),  #  (leopard)
                      ("Pv", "Prionailurusviverrinus_genome", "GCA_018119265.1"),  #  (fishing cat)
                      ("Lq", "Leopardisgeoffroyi_genome", "GCA_018350155.1"),  # Leopardus geoffroyi
                      ("Lc", "Lynxcanadensis_genome", "GCA_007474595.2"),  #  (lynx)
                      ("Aj", "Acinonyxjubatus_genome", "GCA_003709585.1"),  #  (cheetah)
                      ("Pc", "Pumaconcolor_genome", "GCA_003327715.1"),  # (puma)
                      ("Py", "Pumayagouaroundi_genome", "GCA_014898765.1"),  #  (jaguarundi)
                      ("Po", "Pantheraonca_genome", "GCA_004023805.1"),  #  (jaguar)
                      ("Cc", "Caracalcaracal_genome", "GCA_016801355.1"))}  #  (caracal)

    # Carnivora (diverged 50 myo) -> maybe better Caniformes (without cats))
    carnivora_cat_dict = {"Fc": (("Hh", "HyaenaHyaena_genome", "GCA_004023945.1"),  # hyaena
                       ("Ss", "SuricataSuricatta_genome", "GCA_006229205.1"),  # meerkat
                       ("Cg", "CryptoproctaFerox_genome", "GCA_004023885.1"),  # fossa
                       ("Ph", "ParadoxurusHermaphroditus_genome", "GCA_004024585.1"),  # asian palm civet
                       ("Sg", "SpilogaleGracilis_genome", "GCA_004023965.1"),   #  western spotted skunk
                       ("Ll",  "LutraLutra_genome", "GCA_902655055.2"),  # Eurasian river otter
                       ("Pl", "ProcyonLotor_genome",  "GCA_015708975.1"),  # raccoon
                       ("Af", "AilurusFulgens_genome", "GCA_002007465.1"),  # lesser panda
                       ("Zc", "ZalophusCalifornianus_genome", "GCA_009762305.2"),  # sea lion
                       ("Or", "OdobenusRosmarus_genome", "GCA_000321225.1"),  # pacific walrus
                       ("Ml",  "MiroungaLeonina_genome",  "GCA_011800145.1"),  # elephant seal
                       ("Um", "UrsusMaritimus_genome", "GCA_017311325.1"),  # polar bear
                       ("Cf", "CanisFamiliaris_genome", "GCA_000002285.2"))}  # dog

    carnivora_dog_dict = {"Cf": (("Um", "UrsusMaritimus_genome", "GCA_017311325.1"),  # polar bear
                       ("Ml",  "MiroungaLeonina_genome",  "GCA_011800145.1"),  # elephant seal
                       ("Or", "OdobenusRosmarus_genome", "GCA_000321225.1"),  # pacific walrus
                       ("Zc", "ZalophusCalifornianus_genome", "GCA_009762305.2"),  # sea lion
                       ("Af", "AilurusFulgens_genome", "GCA_002007465.1"),  # lesser panda
                       ("Pl", "ProcyonLotor_genome",  "GCA_015708975.1"),  # raccoon
                       ("Ll",  "LutraLutra_genome", "GCA_902655055.2"),  # Eurasian river otter
                       ("Sg", "SpilogaleGracilis_genome", "GCA_004023965.1"),   #  western spotted skunk
                       ("Ph", "ParadoxurusHermaphroditus_genome", "GCA_004024585.1"),  # asian palm civet
                       ("Cg", "CryptoproctaFerox_genome", "GCA_004023885.1"),  # fossa
                       ("Hh", "HyaenaHyaena_genome", "GCA_004023945.1"),  # hyaena
                       ("Ss", "SuricataSuricatta_genome", "GCA_006229205.1"),  # meerkat
                       ("Fc", "CATUS_genome", "GCF_000181335.1"))}  # cat


    # redundant parenthesis allow collecting ref and not ref genomes with the same loop (below)
    mouse_ref_dict = {"Mm": (("Mm", "MUSCULUS_genome", "GCA_000001635.9"))} # GeneBank GRCm39
    rat_ref_dict = {"Rn": (("Rn", "NORVEGICUS_genome", "GCA_000001895.4"))}  # GenBank; Rnor_6.0 -> NOT SAME AS GENOMES
    human_ref_dict = {"Hs": (("Hs", "SAPIENS_genome", "GCA_000001405.28"))} # RefSeq GCF_000001405.39
    cat_ref_dict = {"Fc": (("Fc", "CATUS_genome", "GCF_000181335.1"))}  # pyensembl is using Felis_catu_6.2
    dog_ref_dict = {"Cf": (("Cf", "FAMILIARIS_genome", "GCA_000002285.2"))}  # CanFam3.1

    if ref_species == "Mm" and ref_genome is False:
        genomes_dict = mouse_dict
    elif ref_species == "Mm" and ref_genome is True:
        genomes_dict = mouse_ref_dict
    elif ref_species == "Rn" and ref_genome is False:
        genomes_dict = rat_dict
    elif ref_species == "Rn" and ref_genome is True:
        genomes_dict = rat_ref_dict
    elif ref_species == "Hs" and ref_genome is False:
        genomes_dict = human_dict
    elif ref_species == "Hs" and ref_genome is True:
        genomes_dict = human_ref_dict
    elif ref_species == "Fc" and ref_genome is False:
        genomes_dict = carnivora_cat_dict
    elif ref_species == "Fc" and ref_genome is True:
        genomes_dict = cat_ref_dict
    elif ref_species == "Cf" and ref_genome is False:
        genomes_dict = carnivora_dog_dict
    elif ref_species == "Cf" and ref_genome is True:
        genomes_dict = dog_ref_dict
    else:
        print("Something went wrong")

    # collect genomes
    all_genomes = []
    for ref, species in genomes_dict.items():
        for genome in species:  # e.g. genome = ("Mi", "SPICILEGUS_genome", "GCA_003336285.1")
            all_genomes.append(genome)

    return all_genomes


def map_assembly_contigs(wdir):
    """Makes a dictionary linking genome indexed contig names and pyensembl release contig names - to get input"""

    import pyensembl

    release = 90
    species = "felis catus"

    ensembl = pyensembl.EnsemblRelease(release, species)

    # get contigs available in ensembl release
    contigs = ensembl.contigs()

    # get all available contigs in the NCBI assembly
    with open(wdir + "CATUS_genome_contigs.txt", "r") as f:
        rows = (line.split('\t') for line in f)
        # need index 3 from each list
        all_contigs = {row[0]: row[1:] for row in rows}

    # get contigs from indexed reference assembly
    with open(wdir + "Reference_genomes/CATUS_genome.fasta.fai", "r") as f:
        rows = (line.split('\t') for line in f)
        # need keys
        indexed_contigs = {row[0]: row[1:] for row in rows}

    #mapped_contigs_dict = {}
    #for contig in contigs:  # comes from ensembl release
    #    for indexed_contig in indexed_contigs:  # comes from pybedtools indexed ref genome .fai file
    #        for available_contig, features in all_contigs.items():  # comes from the ncbi assembly file
    #            if available_contig == contig and indexed_contig == features[3]:
    #                mapped_contigs_dict[contig] = indexed_contig
    #            elif contig == features[3]:
    #                mapped_contigs_dict[contig] = features[3]
    #            elif contig == features[1] and "random" not in available_contig:
    #                mapped_contigs_dict[contig] = features[3]

    mapped_contigs_dict = {}
    for contig in contigs:
        for available_contig, features in all_contigs.items():
            if contig == available_contig or contig in features[3]:
                mapped_contigs_dict[contig] = features[5]

    return mapped_contigs_dict



"""




def map_assembly_contigs(wdir):

    # first need to index ref genome (.fai file)

    import pyensembl

    release = 104
    species = "felis catus"

    ensembl = pyensembl.EnsemblRelease(release, species)

    # get contigs available in ensembl release
    contigs = ensembl.contigs()

    # get all available contigs in the NCBI assembly
    with open(wdir + "CATUS_genome_contigs.txt", "r") as f:
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
