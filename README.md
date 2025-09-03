# ![](freeda_logo.png)

FREEDA
======
FREEDA (Finder of Rapidly Evolving Exons in Diverse Assemblies) is a Python built end-to-end 
automated pipeline to detect positive selection, created for experimental cell biologists 
without training in computational biology and molecular evolution. 

- Documentation: [https://ddudka9.github.io/freeda/](https://ddudka9.github.io/freeda/)
- Source code: [https://github.com/DDudka9/freeda](https://github.com/DDudka9/freeda)
- Requirements:
	- min 100 GB of disc space for each taxon e.g. primates (we recommend a hard drive with 500GB space)
	- Stable Internet connection

If you use FREEDA for published work, please cite the original paper:

Dudka D, Akins RB, Lampson MA. FREEDA: An automated computational pipeline guides experimental testing of protein innovation. J Cell Biol. 2023 Sep 4;222(9):e202212084. doi: 10.1083/jcb.202212084. Epub 2023 Jun 26. PMID: 37358475; PMCID: PMC10292211.


FREEDA is published under the GPLv3 license.

Created by Damian Dudka (damiandudka@0gmail.com) and R. Brian Akins

How to get started
==================

MacOS
-----

1. Download the latest MacOS release from the GitHub Releases page: 
	[https://github.com/DDudka9/freeda/releases/tag/v1.0.10-mac](https://github.com/DDudka9/freeda/releases/tag/v1.0.10-mac)
2. Install PyMOL https://pymol.org/2 in your Applications folder (takes minutes)
3. Double-click the .app file to run FREEDA (it might take a minute to load) and follow the steps depicted below.

*Make sure you have Rosetta 2 program installed on MacOS with Apple Silicon chip

(Find exhaustive walkthrough and explanation of the results in Documentation: [https://ddudka9.github.io/freeda/](https://ddudka9.github.io/freeda/))
![](GUI_example.png)

Linux
-----

1. Download the latest Linux release (compatible with Ubuntu 22.04 LTS) from the GitHub Releases page: 
	[https://github.com/DDudka9/freeda/releases/tag/v1.0.10-linux](https://github.com/DDudka9/freeda/releases/tag/v1.0.10-linux)
2. "Right-click" the FREEDA file -> go to Properties -> go to Permissions -> tick Allow executing file as program
3. Press ctrl + alt + t to open Terminal window.
4. Drag and drop the FREEDA file and press ENTER to open FREEDA GUI (it might take a minute to load).

Windows
-------

FREEDA is not available for Windows. However, it is easy to run FREEDA from a free VirtualBox running Linux inside a Windows computer. 
Setting up a VirtualBox takes minutes. See Documentation: [https://ddudka9.github.io/freeda/](https://ddudka9.github.io/freeda/)
