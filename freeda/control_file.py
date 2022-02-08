"""

Makes a control file for PAML analysis

"""


def make_control_file(wdir):
    """ Makes a control file for PAML from the control_file_dict"""

    # make control file for F3x4 model (default)
    with open(wdir + "control_file_F3x4.ctl", "w") as f:
        for index, line in control_file_dict_F3x4.items():
            f.write(line)

    # make control file for F61 model (default)
    with open(wdir + "control_file_F61.ctl", "w") as f:
        for index, line in control_file_dict_F61.items():
            f.write(line)


control_file_dict_F3x4 = \
{0: '      seqfile = input.phy     * sequence data filename\n',
 1: '     treefile = gene.tree   * tree structure file name\n',
 2: '      outfile = output_PAML    * main result file name\n',
 3: '\n',
 4: '        noisy = 0  * 0,1,2,3,9: how much rubbish on the screen\n',
 5: '      verbose = 1  * 0: concise; 1: detailed, 2: too much\n',
 6: '      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic\n',
 7: '                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise\n',
 8: '\n',
 9: '      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs\n',
 10: '    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n',
 11: '*       ndata = 1\n',
 12: '        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n',
 13: '       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\n',
 14: '   aaRatefile =    * only used for aa seqs with model=empirical(_F)\n',
 15: '                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own\n',
 16: '\n',
 17: '        model = 0\n',
 18: '                   * models for codons:\n',
 19: '                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches\n',
 20: '                   * models for AAs or codon-translated AAs:\n',
 21: '                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F\n',
 22: '                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)\n',
 23: '\n',
 24: '      NSsites = 0 1 2 7 8  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;\n',
 25: '                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;\n',
 26: '                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;\n',
 27: '                   * 13:3normal>0\n',
 28: '\n',
 29: '        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below\n',
 30: '        Mgene = 0\n',
 31: '                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n',
 32: '                   * AA: 0:rates, 1:separate\n',
 33: '\n',
 34: '    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated\n',
 35: '        kappa = 2  * initial or fixed kappa\n',
 36: '    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate \n',
 37: '        omega = 1  * initial or fixed omega, for codons or codon-based AAs\n',
 38: '\n',
 39: '    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha\n',
 40: '        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)\n',
 41: '       Malpha = 0  * different alphas for genes\n',
 42: '        ncatG = 8  * # of categories in dG of NSsites models\n',
 43: '\n',
 44: "        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates\n",
 45: ' RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)\n',
 46: '\n',
 47: '   Small_Diff = .5e-6\n',
 48: '    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?\n',
 49: '*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed\n',
 50: '       method = 0  * Optimization method 0: simultaneous; 1: one branch a time\n',
 51: '\n',
 52: '* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,\n',
 53: '* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., \n',
 54: '* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., \n',
 55: '* 10: blepharisma nu.\n',
 56: '* These codes correspond to transl_table 1 to 11 of GENEBANK.\n'}


control_file_dict_F61 = \
{0: '      seqfile = input.phy     * sequence data filename\n',
 1: '     treefile = gene.tree   * tree structure file name\n',
 2: '      outfile = output_PAML    * main result file name\n',
 3: '\n',
 4: '        noisy = 0  * 0,1,2,3,9: how much rubbish on the screen\n',
 5: '      verbose = 1  * 0: concise; 1: detailed, 2: too much\n',
 6: '      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic\n',
 7: '                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise\n',
 8: '\n',
 9: '      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs\n',
 10: '    CodonFreq = 3  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n',
 11: '*       ndata = 1\n',
 12: '        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n',
 13: '       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\n',
 14: '   aaRatefile =    * only used for aa seqs with model=empirical(_F)\n',
 15: '                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own\n',
 16: '\n',
 17: '        model = 0\n',
 18: '                   * models for codons:\n',
 19: '                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches\n',
 20: '                   * models for AAs or codon-translated AAs:\n',
 21: '                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F\n',
 22: '                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)\n',
 23: '\n',
 24: '      NSsites = 0 1 2 7 8  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;\n',
 25: '                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;\n',
 26: '                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;\n',
 27: '                   * 13:3normal>0\n',
 28: '\n',
 29: '        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below\n',
 30: '        Mgene = 0\n',
 31: '                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n',
 32: '                   * AA: 0:rates, 1:separate\n',
 33: '\n',
 34: '    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated\n',
 35: '        kappa = 2  * initial or fixed kappa\n',
 36: '    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate \n',
 37: '        omega = 1  * initial or fixed omega, for codons or codon-based AAs\n',
 38: '\n',
 39: '    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha\n',
 40: '        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)\n',
 41: '       Malpha = 0  * different alphas for genes\n',
 42: '        ncatG = 8  * # of categories in dG of NSsites models\n',
 43: '\n',
 44: "        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates\n",
 45: ' RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)\n',
 46: '\n',
 47: '   Small_Diff = .5e-6\n',
 48: '    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?\n',
 49: '*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed\n',
 50: '       method = 0  * Optimization method 0: simultaneous; 1: one branch a time\n',
 51: '\n',
 52: '* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,\n',
 53: '* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., \n',
 54: '* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., \n',
 55: '* 10: blepharisma nu.\n',
 56: '* These codes correspond to transl_table 1 to 11 of GENEBANK.\n'}