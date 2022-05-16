Installation
============

Mac
---

The easiest way to install FREEDA is to download the prepackaged versions available on the GitHub Releases page. (INSERT LINK) After downloading, double-click the Mac .app file to run FREEDA.


.. _linux installation anchor:

Linux
-----

A prepackaged version of FREEDA is also available for Linux. Download the latest release from the GitHub Releases page. To run FREEDA, it first must be allowed to run as an executable. To do so, right click on FREEDA and select "Properties...". Then, go to the "Permissions" tab and enable "Allow this file to run as a program". Once this is done, double-click the linux .run file to run FREEDA.

Running FREEDA from the Terminal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the above method does not work, FREEDA can also be run from the command line. To do so, open a terminal either by pressing CTRL + ALT + t or by clicking on it in the Applications menu. Then, navigate to the folder where FREEDA was downloaded, allow it to be executed, and run it. For example, if FREEDA was downloaded to the Downloads folder:

.. code-block:: sh

    cd ~/Downloads
    chmod +x freeda
    ./freeda

Replacing ``~/Downloads`` with your download location (note that path names are case sensitive).

After the first time, FREEDA can be run the same way, omitting the ``chmod...`` command.

Windows
-------

FREEDA is not available for Windows. However, it is possible to run FREEDA from a virtual machine running Linux inside a Windows computer. For more information, see :ref:`Setting Up a Virtual Machine <virtual machine anchor>`.
