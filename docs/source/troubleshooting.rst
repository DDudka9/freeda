Troubleshooting
===============

The most common error comes from unstable Internet connection, especially during the first run. Simply re-run the .app.

In case your issue is not covered here please send (**damiandudka0@gmail.com**):
	- print-screen of the GUI
	- *FREEDA-current-date.log* file (in "Raw_data" folder or in the folder indicated by "Set directory" see: :ref:`Usage <Understanding Results>`)
	- (if present) *PAML-current-date.log* file (in "Raw_data" folder or in the folder indicated by "Set directory" see: :ref:`Usage <Understanding Results>`)


Fatal errors (action needed)
-----------------------------------------

**...FATAL ERROR... : Something went wrong (see below). Send screen shot to : Damian Dudka -> damiandudka0@gmail.com**
	
	*REASON:* Follows most unknown errors causing FREEDA to crash.
	
	*SOLUTION:* Please attach a screen shot and "FREEDA-current-date.log" and (if present) "PAML-current-date.log". They will be either in "Raw_data" folder or in the main folder indicated within the GUI ("Set directory").

**...FATAL ERROR... PyMOL not found in the PATH, the Applications folder, or the working directory.**

	*REASON:* This is most likely a iOS user error installing PyMOL outside of the Applications folder. 
	
	*SOLUTION:* Move the PyMOL .app to Applications.

**...FATAL ERROR... : Partial or empty genome file : GENOME.fasta -> likely Internet connection was disrupted -> please run FREEDA again...**

	*REASON:* The most common error at first run - almost always caused by unstable Internet connection or ABORTing the .app when downloading the genomes at first run. 
	
	*SOLUTION:* Simply re-run the app.

**...FATAL ERROR... : Genome : GENOME blast databases failed to build (likely interrupted) -> rerun the pipeline...**

	*REASON:* Most likely the .app crashed or was ABORTed while making local BLAST databases. 
	
	*SOLUTION:* Try re-running the app.

**...FATAL ERROR... : Genome : GENOME blast databases were built partially (likely interrupted) -> rerun the pipeline...**

	*REASON:* Similar to the previous error. FREEDA is expecting a certain size of the generated databases. Partial files may be present when the run was interrupted.
	
	*SOLUTION:* Try re-running the app.

**...FATAL ERROR... : Repetitive coding sequence detected in GENE (min 80bp repeat) -> cannot reliably analyze this GENE**

	*REASON:* FREEDA's ability to find orthologous exons relies on accurate alignment to the reference coding sequence. Repetitive sequences are notoriously difficult to align, therefore FREEDA will not try analyze genes with min. 80bp repeats in coding sequence (e.g. Tacc3 in rodents).
	
	*SOLUTION:* Unfortunately this gene cannot be analyzed reliably. If you still want to try I can run the analysis on my end and let you know if its reliable: damiandudka0@gmail.com

**...FATAL ERROR... : No reliable coding sequence annotation detected for GENE**

	*REASON:* This is usually caused by unreliable annotation of the gene of interest - mostly concerning dog and chicken assemblies.
	
	*SOLUTION:* Make sure there is a known coding sequence for this gene (e.g. visit Ensembl database - `https://useast.ensembl.org/index.html <https://useast.ensembl.org/index.html>`_

**...FATAL ERROR... : Input data generation FAILED for GENE - please remove GENE from analysis -> exiting the pipeline now...**

	*REASON:* Unreliable input generation - usually following another FATAL_ERROR (e.g. no coding transcript found)
	
	*SOLUTION:* Make this is not a pseudogene (e.g. visit Ensembl database - `https://useast.ensembl.org/index.html <https://useast.ensembl.org/index.html>`_
	
**...FATAL ERROR... : At least 3 blast output files contain no matches above threshold : ... for gene name: GENE -> please exclude them and run FREEDA again -> exiting the pipeline now...**

	*REASON:* When using "Common domain expected" option to increase the threshold of BLAST hits (from 60 to 80%), you might end up with no hits for very divergent genes. It might also be due to a gene loss. 
	
	*SOLUTION:* Try the "Exclude species" option using two-letter code for each genome (e.g. Pt Gg for *Pan troglodytes* and *Gorilla gorilla*; see CITATION).


Warnings (no action needed)
---------------------------

**AlphaFold structure not found or not matching Ensembl input collected**

	*REASON:* FREEDA will still run the analysis but without mapping residues onto structure

	.. image:: /_images/GUI_events_No_structure.png

**Questionable alignment of a single exon**

	*REASON:* FREEDA performs additional checks (blue) for a divergent exon and accepts or rejects it
	
	.. image:: /_images/GUI_events_single_exon_warning.png

**Coding sequence is not in frame**

	*REASON:* Either some exons are missing (not the case in example below) or single indels 
	are present (e.g. sequencing errors). FREEDA may either remove this sequence from analysis 
	or remove the indels to force conserved alignment.
	
	.. image:: /_images/GUI_events_CDS_not_in_frame.png

**Early STOP codon detected in final coding sequence**
	
	*REASON:* Coding sequences in the final alignment should not have any STOP codons anymore, 
	except those introduced by imperfect alignment. This is rare and usually occurs in highly 
	divergent genes. FREEDA removes the STOP codon forcing a conserved alignment.
	
	.. image:: /_images/GUI_STOP_codon_warning.png


**Failed check comparing cloned sequence to annotated one for most distant species**
	
	*REASON:* This is a sanity check. Usually <95% identity suggests that alternative exons 
	are used (check supported only for rodents and carnivores).

**Potential misalignment is visible in the protein sequence alignment**

	*REASON:* *Cercopithecus mona* orthologue shows a distinct pattern of non-synonymous 
	substitutions, raising a possibility of misalignment. The user should re-analyze the 
	gene of interest using the "exclude species" option to avoid false positive signature 
	of positive selection.

**Potential misalignment is visible in the raw nucleotide sequence alignment**
	
	*REASON:* *Cercopithecus mona* orthologue shows an unusual out-of-frame 2bp deletion 
	in the middle of the sequence, possibly due to a sequencing error. While FREEDA alignment 
	filtering is robust to prevent a global protein misalignment, a local misalignment can 
	still occur as seen in the protein sequence (see above). The user should re-analyze the 
	gene of interest using the "exclude species" option to avoid false positive signature 
	of positive selection.



Problems with VirtualBox
------------------------

64-bit Operating Systems Won't Show Up
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your Windows system may already be using another Virtual Machine system called HyperV.

Other problems
--------------

PyMOL Tar Unpacking "Cannot Create Symlink"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your data folder may be in a system partition format that does not allow symbolic links. As a work around, you can install PyMol yourself from a package manager, software application, or the open source PyMol repository: `https://github.com/schrodinger/pymol-open-source <https://github.com/schrodinger/pymol-open-source>`.





