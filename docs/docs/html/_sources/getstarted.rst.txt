How to get started
==================

MacOS
------------------

(Compatible with both Intel and Apple Silicon MacOS: Mojave, Big Sur, Monterrey and Ventura)

1. Download the latest MacOS release from the GitHub Releases page: 
	`https://github.com/DDudka9/freeda/releases/tag/v1.0.10-mac <https://github.com/DDudka9/freeda/releases/tag/v1.0.10-mac>`_
2. Install PyMOL `https://pymol.org/2/ <https://pymol.org/2/>`_ in your Applications folder (takes minutes).
3. Double-click the .app file to open FREEDA GUI (it might take a minute to load).

* Make sure you have Rosetta 2 program installed on MacOS with Apple Silicon chip

.. _linux installation anchor:

Linux
-----

(Compatible with Ubuntu 22.04 LTS release `https://ubuntu.com/download/desktop <https://ubuntu.com/download/desktop>`_)

1. Download the latest Linux release from the GitHub Releases page: 
	`https://github.com/DDudka9/freeda/releases/tag/v1.0.10-linux <https://github.com/DDudka9/freeda/releases/tag/v1.0.10-linux>`_
2. "Right-click" the FREEDA file -> go to Properties -> go to Permissions -> tick Allow executing file as program
3. Press ctrl + alt + t to open Terminal window.
4. Drag and drop the extracted file and press ENTER to open FREEDA GUI (it might take a minute to load).


Windows
-----------------------------

FREEDA is not available for Windows. However, it is easy to run FREEDA from a free VirtualBox running Linux inside a Windows computer. Setting up a VirtualBox takes minutes.

.. _virtual machine anchor:

**Setting Up a VirtualBox (10 min)**

One way to set up a virtual machine is:

1. Download a commonly used VirtualBox from `https://virtualbox.org/wiki/Downloads <https://virtualbox.org/wiki/Downloads>`_.

2. Download Linux version Ubuntu 22.04 LTS from `https://ubuntu.com/download/desktop <https://ubuntu.com/download/desktop>`_. In general, operating system files have the suffix ".iso".

3. Open VirtualBox and click on "New" to create a new virtual machine.

   .. image:: /_images/VB1.png

4. Next, enter a name and installation location for your virtual machine. Choose "Linux" for the type and "Ubuntu 22.04 LTS (Jammy Jellyfish) (64-bit)" for the version.

   .. image:: /_images/New_VB2.png

5. On the next screen select disk space. We recommend allocating 500GB to ensure you dont run out of space! YOU CANNOT CHANGE THAT LATER.

   .. image:: /_images/New_VB2_1.png

6. Next allocate RAM and CPU units for the virtual machine. At least 6GB (6144MB) of RAM is recommended in this step (DO NOT allocate the maximum number available!).
   We recommend selecting 2 CPU units to ensure high performance and stability of your Linux (DO NOT allocate the maximum number available!).
   You can change both of those later.

   .. image:: /_images/New_VB3.png

7.  After creating the virtual machine, assign the operating system image to it by clicking the option next to IDE Secondary Device 0 and "Choosing disk file" that you downloaded earlier (.iso file).

   .. image:: /_images/New_VB4.png

   .. image:: /_images/New_VB4_2.png

8. If you selected only 1 CPU earlier, now you can increase it to 2 CPU units by going to Settings -> System -> Processor.

   .. image:: /_images/New_VB4_3.png

9. Open the virtual machine by clicking "Start". Once it has started, choose "Install Ubuntu".

   .. image:: /_images/New_VB5.png

10. Follow the installer instructions. Note that when the screen comes up to "Erase disk and install Ubuntu" this will **not** affect the other files on your computer. It is only erasing the data inside the virtual machine.

   .. image:: /_images/New_VB6.png

11. Continue with the installation. When it is complete, the virtual machine will restart. When it asks to "remove the installation medium", just press ENTER.

   .. image:: /_images/VB7.png

12. If the installation has been completed successfully, the virtual machine should boot to a desktop image.

    .. image:: /_images/New_VB8.png

Finally, we suggest increasing the resolution of the VirtualBox Linux window.
Settings (right top corner) -> Displays -> Resolution 1360x768) to make sure the GUI fits the screen.
We also suggest preventing sleeping -> Settings -> Power -> Screen Blank -> Never to ensure uninterrupted analysis.

Congratulations! To run FREEDA in VirtualBox just follow instructions for (:ref:`Linux`)
-> remember that you will have to download FREEDA within your VirtualBox
(best is to open Mozilla Firefox that is automatically installed by Linux and go to
`https://github.com/DDudka9/freeda/releases/tag/v1.0.10-linux <https://github.com/DDudka9/freeda/releases/tag/v1.0.10-linux>`_)