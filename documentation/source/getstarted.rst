How to get started
==================

MacOS
-----

1. Download the latest MacOS release from the GitHub Releases page: 
	(INSERT LINK)
2. Install PyMOL `https://pymol.org/2/ <https://pymol.org/2/>`_ in your Applications folder
3. Double-click the .app file to open GUI and run FREEDA (it might take a minute to load)


.. _linux installation anchor:

Linux
-----

1. Download the latest Linux release from the GitHub Releases page: 
	(INSERT LINK)
2. Go to the folder and "right-click" the file -> in Properties tick "Allow to run as executable"
3. Double-click the file to open GUI and run FREEDA (it might take a minute to load)


2. Run FREEDA from the command line. To do so, open a terminal either by pressing CTRL + ALT + t or by clicking on it in the Applications menu. 
3. Navigate to the folder where FREEDA was downloaded, allow it to be executed, and run it. 

For example, if FREEDA was downloaded to the Downloads folder copy and enter commands:

.. code-block:: sh

    cd ~/Downloads
    chmod +x freeda_pipeline_GUI
    ./freeda_pipeline_GUI

Otherwise replace ``~/Downloads`` with your download location (path names are case sensitive).

After the first time, FREEDA can be run the same way, omitting the ``chmod...`` command.


Windows
-------

FREEDA is not available for Windows. However, it is easy to run FREEDA from a free Virtual Machine running Linux inside a Windows computer. Setting up a Virtual Machine takes minutes.

.. _virtual machine anchor:

**Setting Up a Virtual Machine (10 min)**

One way to set up a virtual machine is:

1. Download a virtual machine manager. A commonly used program is Oracle's VirtualBox, which can be downloaded from `https://virtualbox.org/wiki/Downloads <https://virtualbox.org/wiki/Downloads>`_.

2. Download a Linux operating system file. There are many options, but a commonly used one is Ubuntu, which can be downloaded from `https://ubuntu.com/download/desktop <https://ubuntu.com/download/desktop>`_. FREEDA requires a Ubuntu version of 16.04.6 or later. In general, operating system files have the suffix ".iso".

3. Open VirtualBox and click on "New" to create a new virtual machine.

   .. image:: /images/VB1.png

4. Enter a name and installation location for your virtual machine. Choose "Linux" for the type and "Ubuntu (64-bit)" for the version. Allocate RAM for the virtual machine. Allocating more RAM will make the virtual machine run faster. At least 4GB of RAM is recommended in this step.

   .. image:: /images/VB2.png

5. On the next screen, create a disk image. 100GB of disk space is needed to run FREEDA on one taxon, but we recommend allocating at least 300GB. The rest of the settings on this screen can remain as their defaults.

   .. image:: /images/VB3.png

6.  After creating the virtual machine, assign the operating system image to it by clicking on the option next to IDE Secondary Device 0 and navigating to the system image file downloaded before.

   .. image:: /images/VB4.png

7. In "Settings" go to "Processor" tab and increase the number of CPUs to 3.



8. Open the virtual machine by clicking "Start". Once it has started, choose "Install Ubuntu".

   .. image:: /images/VB5.png

9. Follow the installer instructions. Note that when the screen comes up to "Erase disk and install Ubuntu" this will **not** affect the other files on your computer. It is only erasing the data inside the virtual machine.

   .. image:: /images/VB6.png

10. Continue with the installation. When it is complete, the virtual machine will restart. When it asks to "remove the installation medium", just press ENTER.

   .. image:: /images/VB7.png

11. If the installation has been completed successfully, the virtual machine should boot to a desktop image.

    .. image:: /images/VB8.png

Congratulations! To run FREEDA now, just follow the steps for a Linux installation at :ref:`Linux Installation <linux installation anchor>`.

12. We suggest increasing the resolution of the VirtualBox Linux window (Settings -> Displays -> Resolution 1360x768) to make sure the GUI fits the screen. We also suggest to prevent sleeping -> Settings -> Power -> Screen Blank -> to ensure robust analysis
