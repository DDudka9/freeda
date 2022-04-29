========
Advanced
========

.. _virtual machine anchor:

Setting Up a Virtual Machine
----------------------------

Although FREEDA is not availabe for Windows, it is possible to run it from a Linux virtual machine set up inside a Windows computer. To set up such a virual machine:

1. Download a virtual machine manager. A commonly used program is Oracle's VirtualBox, which can be downloaded from `https://virtualbox.org/wiki/Downloads <https://virtualbox.org/wiki/Downloads>`.

2. Download a Linux operating system file. There are many options, but a commonly used one is Ubuntu, which can be downloaded from `https://ubuntu.com/download/desktop <https://ubuntu.com/download/desktop>`.

3. Set up the virtual machine. In VirtualBox,

Installing from Source
----------------------

Installing FREEDA from the provided executable as outlined in :doc:`Installation <installation>` is the easiest way to run FREEDA. However, it is also possible to install FREEDA from source files, which allows improved debugging abilities and modification of the source code. To install FREEDA from source:

1. Download the FREEDA source code from GitHub (INSERT LINK). Be sure to download both the main executable, the ``freeda`` package, and the anaconda environment file.

2. Create a conda environment by:

   - Installing Anaconda from a package manager or `https://www.anaconda.com/ <https://www.anaconda.com/>` and following setup instructions.
   - Creating the environment with `conda create --name FREEDA freeda_environment.yml` from a terminal in the directory to which freeda was downloaded.
   - Activating the environment with `conda activate FREEDA`

3. Run FREEDA with ``./freeda_pipeline_GUI``.
