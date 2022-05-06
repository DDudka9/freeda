========
Advanced
========

.. _virtual machine anchor:

Setting Up a Virtual Machine
----------------------------

Although FREEDA is not availabe for Windows, it is possible to run it from a Linux virtual machine set up inside a Windows computer. One way to set up such a virual machine is:

1. Download a virtual machine manager. A commonly used program is Oracle's VirtualBox, which can be downloaded from `https://virtualbox.org/wiki/Downloads <https://virtualbox.org/wiki/Downloads>`_.

2. Download a Linux operating system file. There are many options, but a commonly used one is Ubuntu, which can be downloaded from `https://ubuntu.com/download/desktop <https://ubuntu.com/download/desktop>`_. FREEDA requires a Ubuntu version of 16.04.6 or later. In general, operating system files have the suffix ".iso".

3. Open VirtualBox and click on "New" to create a new virtual machine.

   .. image:: //home/bakins/freeda/docs/images/VB1.png

4. Enter a name and installation location for your virtual machine. Choose "Linux" for the type and "Ubuntu (64-bit)" for the version. Allocate RAM for the virtual machine. Allocating more RAM will make the virtual machine run faster. At least 4GB of RAM is recommended in this step.

   .. image:: //home/bakins/freeda/docs/images/VB2.png

5. On the next screen, create a disk image. At least 200GB of disk space is required to store FREEDA and its genomes, but 500GB is recommended. The rest of the settings on this screen can remain as their defaults.

   .. image:: //home/bakins/freeda/docs/images/VB3.png

6.  After creating the virtual machine, assign the operating system image to it by clicking on the option next to IDE Secondary Device 0 and navigating to the system image file downloaded before.

   .. image:: //home/bakins/freeda/docs/images/VB4.png

7. Open the virtual machine by clicking "Run". Once it has started, choose "Install Ubuntu".

   .. image:: //home/bakins/freeda/docs/images/VB5.png

8. Follow the installer instructions. Note that when the screen comes up to "Erase disk and install Ubuntu" this will **not** affect the other files on your computer. It is only erasing the data inside the virtual machine.

   .. image:: //home/bakins/freeda/docs/images/VB6.png

9. Continue with the installation. When it is complete, the virtual machine will restart. When it asks to "remove the installation medium", just press ENTER.

   .. image:: //home/bakins/freeda/docs/images/VB7.png

10. If the installation has been completed successfully, the virtual machine should boot to a desktop image.

    .. image:: //home/bakins/freeda/docs/images/VB8.png

Congratulations! To run FREEDA now, just follow the steps for a Linux installation at :ref:`Linux Installation <linux installation anchor>`.

   

Installing from Source
----------------------

Installing FREEDA from the provided executable as outlined in :doc:`Installation <installation>` is the easiest way to run FREEDA. However, it is also possible to install FREEDA from source files, which allows improved debugging abilities and modification of the source code. To install FREEDA from source:

1. Download the FREEDA source code from GitHub (INSERT LINK). Be sure to download both the main executable, the ``freeda`` package, and the anaconda environment file.

2. Create a conda environment by:

   - Installing Anaconda from a package manager or `https://www.anaconda.com/ <https://www.anaconda.com/>`_ and following setup instructions.
   - Creating the environment with `conda create --name FREEDA freeda_environment.yml` from a terminal in the directory to which freeda was downloaded.
   - Activating the environment with `conda activate FREEDA`

3. Run FREEDA with ``./freeda_pipeline_GUI``.
