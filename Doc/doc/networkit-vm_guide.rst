===============================
NetworKit VM Installation Guide
===============================

This document will guide you through the installation of NetworKit in a virtual machine.

.. note:: This type of installation is only recommended for users running Microsoft Windows or testing NetworKit. The
  performance is bounded by the virtual machine.


Step 1 - Installing Oracle VM VirtualBox
----------------------------------------

If you do not already have VirtualBox installed on your system, visit https://www.virtualbox.org/wiki/Downloads and download the right VirtualBox
binary for your system. Once downloaded, run the setup and follow the instructions to install VirtualBox on your system.


Step 2 - Import NetworKit VM into VirtualBox
--------------------------------------------

Download the NetworKit VM from TODO. After that, open VirtualBox and click on „File -> Import Appliance ...“ or the
corresponding entry in your system language.

.. image:: resources/networkit_vm_import.png
  :width: 550
  :align: center

|

In the opening dialog, click on the small folder icon and specify the path to the NetworKit VM file you downloaded.
Hitting „Next“ will show you some options for setting up the virtual machine (number of cores, amount of RAM, etc.).
The standard settings should be fine for the moment and you can finish this step by clicking the „Import“ button.


Step 3 - Running the NetworKit VM
---------------------------------

After the previous steps, you will find an entry „NetworKit“ in the left list of available virtual machines in VirtualBox.
To start the virtual machine, double click on „NetworKit“ or select the entry and click on the „Start“ button in the top
menu bar.

.. image:: resources/networkit_vm_start.png
	:align: center


|

After booting up the virtual machine, you will be prompted to enter a password which is „networkit“. The current version of
NetworKit is preinstalled and can be found in the folder linked on the desktop.
