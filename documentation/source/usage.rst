Usage
=====

Using FREEDA GUI
----------------

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
5. *...WARNING:...* - suboptimal behavior - no action needed (see :ref:`Troubleshooting`)
6. *...FATAL_ERROR:...* - critical failure - action needed (see :ref:`Troubleshooting`)
7. Putting your computer to sleep should not interfere with the analysis (will resume after awaking)
8. We recommend the free software UNIPRO by Ugene for viewing alignment files `http://ugene.net/download-all.html <http://ugene.net/download-all.html>`_
9. When opening structure files (.pse) with PyMOL you can click "Skip activation" - no license is ever needed to view the structure files 


Understanding Events Window
---------------------------

**FIRST RUN ONLY:**

1. If running on MacOS, FREEDA will prompt the user to download PyMOL from `https://pymol.org/ <https://pymol.org/>`_
2. FREEDA will first download Ensembl release information for a chosen taxon needed for input extraction (lots of text)
3. Next FREEDA will download and decompress the reference genome for selected taxon (this will take several minutes depending on Internet speed)
4. Next FREEDA will download, decompress genomic related assemblies and create local blast databases for each (allow 1-2h for this step dependent on Internet speed)


	.. image:: /images/GUI_genomes_downloaded.png
   
   There will be 19-21 genomic assemblies downloaded for each selected taxon.
   This step is triggered each time you select a new working directory ("Set directory" in GUI)!

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
	
**Exon calling - this contig likely contains a retro-duplication**

	.. image:: /images/GUI_events_RETRO.png

**Exon calling - this contig is missing the last two syntenic exons expected**

	.. image:: /images/GUI_events_syntenic_5prime.png

**Exon calling - this contig contains only the last two syntenic exons expected**

	*is MISSING* and *does not have intron* are functionally equivalent - syntenic exon not found

	.. image:: /images/GUI_events_syntenic_3prime.png

**Resolution of very recent duplications (or heterozygous loci)**
	
	*This step is triggered only when at least two contigs bear the same number of likely syntenic exons (e.g. very recent segmental duplications). If the likelihood of synteny is the same - each exon will be compared to the corresponding reference exon using a hamming distance algorithm. The contig with the lowest hamming distance is selected as the most likely orthologous locus (most conserved)*
	
	(RERUN MICB in primates)
		
	"RETRO_score" is always active and flags retro-duplications.
	
	"Synteny_score" is enabled only when "Duplication expected" advanced option is selected within the GUI.
	
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

	*Gene* - Gene name provided
	
	*Pos. select.* - Is there evidence of positive selection acting on the gene?
	
	*LRT* - Likelihood Ratio Test that determines statistical support for evidence of positive selection (>5.99 -> p=0.05)
	
	*p-value* - Directly linked to the LRT value
	
	*CDS cover.* - Percentage of codons analyzed as compared to the reference coding sequence (microexons are excluded from this calculation)
	
	*species* - Number of species (orthologues) analyzed. **Less than 6 species may yield unreliable results**
	
	*pr < 0.9* - Number of all residues that might be evolving under positive selection
	
	*pr >= 0.9* - Number of residues with high probability of being under positive selection
	

**Folder with all results (inside user indicated "Set directory")**

	.. image:: /images/Working_directory_Raw_data.png
	.. image:: /images/Working_directory_Results_data.png

**Exemplary nucleotide alignment (open in UNIPRO Ugene)**

	*Cenpo_raw_nucleotide_alignment.fasta*

	.. image:: /images/Exemplary_nucleotide_alignment.png
	
	Marked is an indel (likely deletion in *Apodemus sylvaticus*) before any processing. Region marked will be removed as it cannot be analyzed. Inspect this file to find which species causes loss of regions from final alignment. 

**Exemplary protein alignment (opne in UNIPRO Ugene)**

	*Cenpo_protein_alignment.fasta*

	.. image:: /images/Exemplary_protein_alignment.png
	
	Marked is the same indel (see above) after it has been processed. Although only 9bp are missing, they span 4 codons. Therefore 4 amino acids were removed from each species (including the first species - after the analysis is complete, FREEDA adds back the missing amino acids to show what was removed). Inspect this file for frameshifts. Use abbreviations displayed here to exclude species.

**Exemplary gene tree (open in Figtree)**
	
	*Cenpo.tree*
	
	.. image:: /images/Exemplary_gene_tree.png

**Results worksheet**

	*PAML_result-10-31-2022-13-02_F3X4*

	.. image:: /images/Exemplary_results_sheet.png
	
	Here you can find probability of positive selection acting on each recurrently changing residue (displayed on top).
	
**Residues under positive selection mapped onto referene CDS**
	
	*Cenpo_PAML_graph_F3X4.tif*
	
	.. image:: /images/Exemplary_graph.png
	
	Top graph (black) shows recurrently changing residues. Middle graph (blue) shows residues that evolve under positive selection with more or less probability (0-7-1.0). Bottom graph (magenta) shows residues with the highest probability of evolving under positive selection. Gray regions have been excluded from analysis (e.g. indels).

**Residues under positive selection mapped onto structural prediction (open in PyMOL)**

	*Cenpo_Mm.pse*

	.. image:: /images/Exemplary_protein_structure.png
	
	You can rotate the structure to have a better look at the position of each residue under positive selection. For details on how to further analyze your structure see: PyMOL wiki `https://pymolwiki.org/index.php/Practical_Pymol_for_Beginners <https://pymolwiki.org/index.php/Practical_Pymol_for_Beginners>`_ and useful user guide `https://pymol.sourceforge.net/newman/userman.pdf <https://pymol.sourceforge.net/newman/userman.pdf>`_
	
	



