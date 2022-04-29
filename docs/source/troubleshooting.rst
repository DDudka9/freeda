Troubleshooting
===============

Problems with FREEDA
--------------------

File Not Found during Genome Downloading and Unpacking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A genome download may fail for a number of reasons, including an interruption in the internet connection. To fix this, navigate to the "Genomes" folder in your chosen working directory and **delete the empty genome file** (it will likely be a .fasta file with a size of 0 bytes). Rerunning FREEDA after this will usually resolve the issue.

PyMol Tar Unpacking "Cannot Create Symlink"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your data folder may be in a system partition format that does not allow symbolic links. As a work around, you can install PyMol yourself from a package manager, software application, or the open source PyMol repository: `https://github.com/schrodinger/pymol-open-source <https://github.com/schrodinger/pymol-open-source>`.


Problems with VirtualBox
------------------------

64-bit Operating Systems Won't Show Up
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your Windows system may already be using another Virtual Machine system called HyperV.
