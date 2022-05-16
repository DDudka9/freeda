Usage
=====

Using FREEDA from the GUI
-------------------------

**INSERT GUI IMAGE**

FREEDA's graphical interface allows users to analyze up to five proteins at a time with a choice of evolutionary clade and codon frequency model. To use FREEDA:

1. Choose from rodents, primates, carnivores, or pheasants to analyze.

2. Choose a codon frequency model. F3X4 is the default model, while F61 is more conservative. If both models are used, only sites that are identified by both models are reported.

3. Enter up to five gene names to analyze. These gene names should match the UniProt entries for the protein and clade chosen. For example, Centromere Protein C would be entered as Cenpc for rodents and as CENPC for primates. If you are unsure of the gene name for your protein of interest, search on the `UniProt Website <https://www.uniprot.org>`_ for your protein of interest and copy the section labeled "Gene".

4. Set the working directory. This is the folder where FREEDA will download genomes, analyze them, and save the results. This directory will be large because of the downloaded genomes, so choose a location with at least 200GB of free space.

5. Click "Analyze".

6. FREEDA will start to run and the activity log will show up in the Events window on the right side of interface. The first time FREEDA is run on a given clade it will require a few hours to download genomes. This step only needs to be performed once, though, and subsequent runs will reuse the downloaded genomes. After genome downloading, FREEDA's analysis takes an hour or less per protein for most proteins. Proteins with extremely large genes, such as those with large introns, may take longer.

7. Once FREEDA is finished, an overview of the results can be found in the Results window. This window will show whether each gene shows signals of positive selection and data about this conclusion. The analyzed data, results, and structures can be found in the working directory chosen in step 4.

Optional Controls
^^^^^^^^^^^^^^^^^

FREEDA offers additional options to control the analysis. They can be used to include already-known information when running FREEDA or to fine tune analysis from previous FREEDA results.

1. If there is a duplication expected of one of the genes, check the "Duplication expected" box, which will <WHAT EXACTLY DOES THIS DO>.

2. If there are known residues or regions in one of the genes, they can be labeled. The "start" and "end" columns mark amino acid residue numbers (where 1 is the N-terminal residue). The "label" column can be any string. These labels do not affect the FREEDA analysis, but they will be shown on the final PyMOL structure.

3. If there are problems with specific genomes, they can be excluded with the "Exclude species" box. The input should be two letters for each species with spaces in between (for example, "Mp Rr" would exclude Mus pahari and Rattus rattus).

Using FREEDA from the Command Line
----------------------------------


