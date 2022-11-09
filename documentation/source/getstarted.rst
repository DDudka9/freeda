How to get started
==================

MacOS
-----

1. Download the latest MacOS release from the GitHub Releases page: 
	(INSERT LINK)
2. Install PyMOL `https://pymol.org/2/ <https://pymol.org/2/>`_ in your Applications folder (takes minutes)
3. Double-click the .app file to open FREEDA GUI (it might take a minute to load)


.. _linux installation anchor:

Linux
-----

(tested on Ubuntu 22.04.1 LTS release `https://ubuntu.com/download/desktop <https://ubuntu.com/download/desktop>`_)

1. Download the latest Linux release from the GitHub Releases page: 
	(INSERT LINK)
2. "Right-click" the file -> "Properties" -> "Permissions" -> tick "Allow executing file as program"
3. Press "ctrl"+ "alt" + "t" to open Terminal window. Drag and drop FREEDA app file and press ENTER to open FREEDA GUI (it might take a minute to load)


Windows
-------

FREEDA is not available for Windows. However, it is easy to run FREEDA from a free Virtual Machine running Linux inside a Windows computer. Setting up a Virtual Machine takes minutes.

.. _virtual machine anchor:

**Setting Up a Virtual Machine (10 min)**

One way to set up a virtual machine is:

1. Download a virtual machine manager. A commonly used program is Oracle's VirtualBox, which can be downloaded from `https://virtualbox.org/wiki/Downloads <https://virtualbox.org/wiki/Downloads>`_.

2. Download Linux version Ubuntu 22.04.1 LTS from `https://ubuntu.com/download/desktop <https://ubuntu.com/download/desktop>`_. In general, operating system files have the suffix ".iso".

3. Open VirtualBox and click on "New" to create a new virtual machine.

   .. image:: /images/VB1.png

4. Enter a name and installation location for your virtual machine. Choose "Linux" for the type and "Ubuntu (64-bit)" for the version. Allocate RAM for the virtual machine. Allocating more RAM will make the virtual machine run faster. At least 4GB (4096MB) of RAM is recommended in this step.

   .. image:: /images/VB2.png

5. On the next screen, create a disk image. 100GB of disk space is needed to run FREEDA on one taxon, but we recommend allocating 500GB to be able to run all taxons. The rest of the settings on this screen can remain as their defaults.

   .. image:: /images/VB3.png

6.  After creating the virtual machine, assign the operating system image to it by clicking on the option next to IDE Secondary Device 0 and choosing the disk file you downloaded before (.iso file).

   .. image:: /images/VB4.png

7. We recommend increasing performance by allocating at least 2 CPU units by going to "Settings" -> "System" -> "Processor".


8. Open the virtual machine by clicking "Start". Once it has started, choose "Install Ubuntu".

   .. image:: /images/VB5.png

9. Follow the installer instructions. Note that when the screen comes up to "Erase disk and install Ubuntu" this will **not** affect the other files on your computer. It is only erasing the data inside the virtual machine.

   .. image:: /images/VB6.png

10. Continue with the installation. When it is complete, the virtual machine will restart. When it asks to "remove the installation medium", just press ENTER.

   .. image:: /images/VB7.png

11. If the installation has been completed successfully, the virtual machine should boot to a desktop image.

    .. image:: /images/VB8.png

Congratulations! To run FREEDA now, just follow the steps for a Linux installation at :ref:`Linux Installation <linux installation anchor>`.

12. We suggest increasing the resolution of the VirtualBox Linux window (Settings -> Displays -> Resolution 1360x768) to make sure the GUI fits the screen. We also suggest to prevent sleeping -> "Settings" -> "Power" -> "Screen Blank" -> "Never" to ensure uninterrupted analysis.

To run FREEDA in Virtual Machine just follow instructions for (ref:`Linux`)
