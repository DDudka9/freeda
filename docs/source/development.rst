Development
===========

Pyinstaller
-----------

FREEDA is packaged using PyInstaller (`documentation <https://pyinstaller.org/>`_). PyInstaller automatically finds the Python interpreter and relevant packages in a system and freezes them for distribution. However, PyInstaller does not find all the required packages for FREEDA automatically. In general, packages which are available only through Conda, and not PyPI, will not be automatically packaged by PyInstaller.

Manually Adding PyInstaller Packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To ensure these packages are included, their binaries and associated library files must be moved to the "include" folder. To find which files are needed:

    1. Build a new conda environment containing only the package(s) to add using ``conda create --name include <LIST PACKAGE NAMES>``, replacing ``<LIST PACKAGE NAMES>`` with any desired packages.

    2. Find the environment folder generated at "<YOUR_CONDA>/envs/include/"

    3. Copy any used binary files from the "bin/" environment folder to "freeda/include/bin/" in the development folder.

    4. Copy all library files from the "lib/" environment folder to "freeda/include/lib/" in the devevelopment folder.

    5. Check that any subprocess calls in the FREEDA source code are surrounded by ``pyinstaller_compatibility.resource_path(<COMMAND>)``.

Including tkinter
^^^^^^^^^^^^^^^^^

TKinter is automatically installed on MacOS systems. For Linux systems, you may need to manually install tkinter with ``sudo apt-get install python3-tk`` (or the equivalent for your package manager and python version).

Setting up the Virtual Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PyInstaller works best when it is run from a virtual environment containing only the necessary files. Running PyInstaller from a conda environment will create much larger output files. Create a virtual enviroment folder and activate it with::
    
    python3 -m venv venv
    source venv/bin/activate

Note that the activate script must be run with "source" before it.

Install the required packages with pip using::

    pip install matplotlib numpy scipy pybedtools 

If the build fails for pybedtools, download the lastest python devolepment tools from your package manager with::

    sudo <PACKAGE-MANAGER> install python3.8-dev

Replacing ``<PACKAGE-MANAGER>`` with your chosen one and replacing ``python3.8-dev`` with the version of python you are using.

Installing UPX
^^^^^^^^^^^^^^

PyInstaller can reduce the size of compiled executables using `UPX (the Ultimate Packer for eXecutables) <https://upx.github.io/>`_. Using UPX is easy: just download it (either from the link or through your package manager), then run PyInstaller as normal.


Creating a .spec File 
^^^^^^^^^^^^^^^^^^^^^
PyInstaller uses files with the .spec extension. These are Python files used to tell the program what to build and which files to include. These files should already be included at one-


Documentation
-------------

This FREEDA documentation is written using `Sphinx <https://sphinx-doc.org/>`_. Sphinx generates HTML documentation from multiple files like this one. Files are related using Table of Contents trees, like the one found in the "index.rst" file. New documentation files can be added by creating a new file with the ".rst" extension in the "source" folder of the docs. The name of this file can then be added to a Table of Contents for access. More information can be found online at the `Sphinx tutorial
<https://www.sphinx-doc.org/en/master/tutorial/index.html>`_.

Sphinx documentation files can be written either in `Markdown <https://www.markdownguide.org/>`_, like the README file, or in `ReStructuredText <https://www.writethedocs.org/guide/writing/reStructuredText/>`_, like this file.

