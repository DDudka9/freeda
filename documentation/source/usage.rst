Usage
=====

Using FREEDA GUI
-------------------------

	.. image:: /images/GUI_input.png

FREEDA's graphical user interface allows analyzing up to 5 genes at a time.

**Mandatory steps:**

a. Provide gene name (when in doubt consult: `https://www.uniprot.org/ <https://www.uniprot.org/>`_)
b. Select reference species
c. Click "Set directory" (we recommend this to be a hard drive)
d. Click "Analyze"


**Optional steps:**

e. Select *Advanced Options* (default OFF):
	
	(We recommend visiting Ensembl or Uniprot databases before using)
	
	- *Duplication expected* (increases confidence in finding true orthologues but decreases recovery of full coding sequence by forcing that each exon must be syntenic at both ends)

	- *Tandem duplication expected* (essential option for tandem duplications but decreases recovery of full coding sequence by restricting exon search to regions closest to exact blast hits)
		
	- *Long introns expected* (increases recovery of full coding sequence but slows down the analysis by searching for exons further away from exact blast hits)
		
	- *Common domains expected* (speeds up the analysis but may miss divergent exons by raising similarity threshold for blast from 60% to 80%)
		
f. *Label residues* (the region of interest will be labeled and colored in the protein structure)

g. *Codon frequencies* (default F3X4)
	
	(We recommend consulting PAML User Guide: `http://abacus.gene.ucl.ac.uk/software/paml.html <http://abacus.gene.ucl.ac.uk/software/paml.html>`_)
	
h. *Exclude species* (recommended in case of observing frameshifts in the alignment; use abbreviation next to a species causing frameshit e.g. Gd for *Grammomys dolichurus*)

i. *ABORT* (generally not recommended but useful to instantly stop the analysis - it shuts down the whole app)



**TIPS:**

1. Stable Internet connection allows smooth download of all genomes (see :ref:`Troubleshooting`)
2. Make a dedicated folder for FREEDA (and always use the same one to avoid downloading genomes again)
3. You can scroll text in the Events Window up-down and left-right
4. *...NOTE:...* - information for the user - no action needed
5. *...WARNING:...* - suboptimal behavior (e.g. no matching structure prediction found) - no action needed (see :ref:`Troubleshooting`)
6. *...FATAL_ERROR:...* - critical failure (e.g. no Internet connection) - FREEDA will stop (see :ref:`Troubleshooting`)
7. Putting your computer to sleep should not interfere with the analysis (will resume after awaking)
8. We recommend the free software UNIPRO by Ugene for viewing alignment files `http://ugene.net/download-all.html <http://ugene.net/download-all.html>`_
9. When opening structure files (.pse) with PyMOL you can click "Skip activation" - no license is ever needed to view the structure files 


Understanding Events Window
-------------------------------

**FIRST RUN ONLY:**

1. If running on MacOS, FREEDA will prompt the user to download PyMOL from `https://pymol.org/ <https://pymol.org/>`_
2. FREEDA will first download Ensembl release information for a chosen taxon needed for input extraction (lots of text)
3. Next FREEDA will download and decompress the reference genome for selected taxon (this will take several minutes depending on Internet speed)
4. Next FREEDA will download, decompress genomic assemblies relevant for selected taxon and create local blast databases for each (allow 1-2h for this step dependent on Internet speed)


	.. image:: /images/GUI_genomes_downloaded.png
   
   There will be 19-21 genomic assemblies downloaded for each selected taxon.
   This step is triggered each time you select a new working directory!

**EACH TIME**

**Checking for genomic assemblies and input extraction**

	.. image:: /images/GUI_input_extraction.png
	
**Searching for homologous sequences using blast**

	.. image:: /images/GUI_events_tblastn.png

**Parsing blast results into separate contigs**

	.. image:: /images/GUI_events_analyzing_gene.png

**Initial alignment of each contig**

	.. image:: /images/GUI_events_aligning_contigs.png
	
**Exon calling - this contig does not contain any syntenic exons expected**

	.. image:: /images/GUI_events_no_introns.png
	
**Exon calling - this contig contains all syntenic exons expected**

	.. image:: /images/GUI_events_syntenic_5_3.png
	
**Exon calling - this contig likely contains a retro-copy of the coding sequence**

	.. image:: /images/GUI_events_RETRO.png

**Exon calling - this contig is missing the last two syntenic exons expected**

	.. image:: /images/GUI_events_syntenic_5prime.png

**Exon calling - this contig contains the last two syntenic exons expected**

	*is MISSING* and *does not have intron* are functionally equivalent - syntenic exon not found

	.. image:: /images/GUI_events_syntenic_3prime.png
	
**Validating single syntenic exons cloned from selected contigs**

	Additional checks are performed if alignment score <0.75; exon is rejected if alignment score <0.60

	.. image:: /images/GUI_events_single_exons.png

**Detecting positive selection**
	
	*Analysis completed* - time it took to find orthologous exons for all analyzed genes. 
	Final multiple sequence alignment is then made for the first gene. Coding sequences with 
	Alignment score <0.69 are eliminated as either containing frameshifts or missing too many exons. 
	Phylogenetic tree for the gene is made based on the nucleotide alignment. 
	PAML analysis starts for the first gene.
	
	.. image:: /images/GUI_events_Analysis_completed.png

	.. image:: /images/GUI_events_LRTs.png
	


Understanding Results
---------------------

**Quick look up table within the GUI**

	.. image:: /images/GUI_result_table.png

**Working directory folder**

**Exemplary nucleotide alignment**

**Exemplary protein alignment**

**Results worksheet**

**Residues under positive selection mapped onto referene CDS**

**Residues under positive selection mapped onto structural prediction**




Troubleshooting
---------------

In case your issue is not covered here please send print-screen 
and *FREEDA-current-date.log* or *PAML-current-date.log* files ("Raw_data" folder)
to **damiandudka0@gmail.com**

**Unstable Internet connection when downloading genomes**
	
	*SOLUTION*: Secure Internet connection, close app, run again (finished downloads will not be affected)
	
	.. image:: /images/GUI_genome_download_error.png

**AlphaFold structure not found or not matching Ensembl input collected**

	*SOLUTION*: No action needed, FREEDA will still run the analysis (without structure output)

	.. image:: /images/GUI_events_No_structure.png

**Questionable alignment of a single exon**

	*SOLUTION*: No action needed, FREEDA performs additional checks (blue) and accepts or rejects the exon
	
	.. image:: /images/GUI_events_single_exon_warning.png

**Coding sequence is not in frame**

	*SOLUTION*: No action needed, likely either some exons are missing (not the case in example below) 
	or single indels (e.g. sequencing errors) - FREEDA will likely remove this sequence from analysis 
	or remove the indels to force conserved alignment
	
	.. image:: /images/GUI_events_CDS_not_in_frame.png

**Failed check comparing cloned sequence to annotated one for most distant species (only rodents and carnivores)**
	
	*SOLUTION*: No action needed, this is a sanity check - usually <95% identity suggests alternative exons used
	



** **
