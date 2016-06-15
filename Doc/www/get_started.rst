.. |separator| raw:: html

	<div style="padding-top: 25px; border-bottom: 1px solid #d4d7d9;"></div>


.. _get_started:

===========
Get Started
===========

We support three ways to install NetworKit:

- `NetworKit Virtual Machine`_: Download and try NetworKit preinstalled on a virtual machine. This is strongly recommended for users using Microsoft Windows.

- `Pip install`_: Download the NetworKit Python package with pip. This is the easier way to get NetworKit but you can only use NetworKit via Python this way.

- `Build NetworKit from Source`_: Clone or download the source code of NetworKit and build the C++ and Python modules from source.



With NetworKit as a Python extension module, you get access to native high-performance code and can at the same time work interactively in the Python ecosystem. Although the standard Python interpreter works fine, we recommend `IPython <http://ipython.readthedocs.org/en/stable/>`_ as a great environment for scientific workflows. View the `IPython Quickstart Guide`_ for installation instructions and how to use NetworKit with IPython.


Once you have installed NetworKit, please make sure to check out our `NetworKit UserGuide <http://nbviewer.ipython.org/urls/networkit.iti.kit.edu/data/uploads/docs/NetworKit_UserGuide.ipynb>`_ for an overview of the features provided in NetworKit.

|separator|

.. _NetworKit Virtual Machine:

Install the NetworKit Virtual Machine
=====================================

If you want a quick and easy way to try NetworKit for your purposes or you use a Microsoft Windows operating system, we strongly recommend the installation of our NetworKit virtual machine.

A detailed installation guide can be found `here <_static/Installation-Guide.pdf>`_.


|separator|


.. _Pip install:

Install NetworKit via Pip
=========================

.. _Python Requirements:

Requirements
~~~~~~~~~~~~

You will need the following software to install NetworKit as a python package:

- A modern C++ compiler, e.g.: `g++ <https://gcc.gnu.org>`_ (>= 4.8) or `clang++ <http://clang.llvm.org>`_ (>= 3.7)
- Python 3 (>= 3.4 is recommended, 3.3 supported)
- `Pip <https://pypi.python.org/pypi/pip>`_
- `SCons <http://scons.org>`_: Please note that SCons is only available for Python 2. For installation via pip, we have a script that builds the C++ part of NetworKit, so you can try it without SCons.

NetworKit uses some additional external Python packages. While you do not need them to run NetworKit, it is strongly recommended to install them in order to use all the features of NetworKit:

- scipy
- numpy
- readline
- matplotlib
- networkx
- tabulate

You can use the command :code:`pip3 install scipy numpy readline matplotlib networkx tabulate` on your terminal to install all packages at once. During the installation of NetworKit, the setup will check if the external packages NetworKit uses are available and print warnings at the end of the installation process. If you do not see any warnings, your system should be ready to use NetworKit.


Install NetworKit
~~~~~~~~~~~~~~~~~

Run :code:`[sudo] pip[3] install [--user] networkit` from your command line to install the Python package *networkit*.

You can remove NetworKit completely by using the command :code:`[sudo] pip[3] uninstall networkit`.


|separator|

.. _Build NetworKit from Source:

Build NetworKit from Source
===========================

You can clone NetworKit from `AlgoHub <http://algohub.iti.kit.edu/parco/NetworKit/NetworKit/>`_ with Mercurial or download the source code as a `Zip file <https://networkit.iti.kit.edu/data/uploads/networkit.zip>`_.

Requirements
~~~~~~~~~~~~

You will need the following software to install NetworKit as a Python package:

- A modern C++ compiler, e.g.: `g++ <https://gcc.gnu.org>`_ (>= 4.8) or `clang++ <http://clang.llvm.org>`_ (>= 3.7)
- `SCons <http://scons.org>`_: Please note that SCons is only available for Python 2. For the different build targets, SCons is mandatory.
- `Google Test <https://github.com/google/googletest>`_ (only needed if you want to build the unit tests, which is recommended)

Building NetworKit
~~~~~~~~~~~~~~~~~~

This section describes how to build NetworKit including the Python functionality. If you do not wish to install NetworKit as a Python package, please refer to `Building Only the C++ Core`_.

For building NetworKit including the Python functionality, make sure to also install the software from the `Python Requirements`_ listed in the `Pip install`_.

After all requirements are installed, switch to the top folder of NetworKit and run the script *setup.py* with the following options:
::
	python3 setup.py build_ext --inplace [--optimize=V] [-jX]

The script will call SCons to compile NetworKit as a library and then build the extensions in the folder *src/python*. By default, NetworKit will be built with the amount of available cores in optimized mode. It is possible to add the options :code:`--optimize=V` and :code:`-jN` the same way it can be done to a manual SCons call, to specify the optimization level and the number of threads used for compilation. The setup script provides more functionality and can be used with pip aswell:
::
	pip3 install -e ./

will compile NetworKit, build the extensions and on top of that temporarily install NetworKit so that it is available on the whole system. This can be undone by calling :code:`pip3 uninstall networkit`.
::
	python3 setup.py clean [--optimize=V]

will remove the extensions and its build folder as well as call SCons to remove the NetworKit library and its build folder specified by :code:`--optimize=V`.

Note: All of the above installation command may require root privileges depending on your system, so try this accordingly. If you do not have root privileges, add :code:`--user` to your command.


Building Only the C++ Core
~~~~~~~~~~~~~~~~~~~~~~~~~~

In case you do not need NetworKit's Python functionality, this section describes how to build the C++ parts only.

We recommend SCons for building the C++ part of NetworKit. Individual settings for your environment will be read from a configuration file. As an example, the file *build.conf.example* is provided. Copy this to *build.conf* and edit your environment settings. Then call Scons.

The call to SCons has the following options:
::
	scons --optimize=<level> --target=<target>

where :code:`<level>` can be

- :code:`Dbg` debug
- :code:`Opt` optimized
- :code:`Pro` profiling

and :code:`target` can be

- :code:`Core` build NetworKit as a library, required for the Python extenstion through Cython.
- :code:`Tests` build executable for the unit tests (requires GoogleTest).
- :code:`Lib` build NetworKit as a library and create symbolic links.

For example, to build NetworKit as an optimized library, run
::
	scons --optimize=Opt --target=Lib

To speed up the compilation on a multicore machine, you can append :code:`-jX` where *X* denotes the number of threads to compile with.

Logging is enabled by default. If you want to disable logging functionality, add the following to your scons call:
::
	--logging=no

Test
~~~~

You actually do not need to build and run our unit tests. However, if you experience any issues with NetworKit, you might want to check, if NetworKit runs properly. Please refer to the :ref:`devGuide-unitTests` section in our :ref:`devGuide`.


|separator|


.. _IPython Quickstart Guide:

Use NetworKit with IPython
==========================

First make sure you have installed IPython, e.g. via pip: :code:`pip3 install ipython`.

IPython Terminal
~~~~~~~~~~~~~~~~

If you want to use NetworKit in the IPython terminal, type the following commands in your OS terminal:
::
	ipython3

.. code-block:: python

	from networkit import *

The first line opens the IPython terminal. The second line imports the *networkit* Python module. After that, you should be able to use NetworKit interactively. For usage examples, refert to the `NetworKit UserGuide <http://nbviewer.ipython.org/urls/networkit.iti.kit.edu/data/uploads/docs/NetworKit_UserGuide.ipynb>`_.

IPython Notebook/jupyterhub
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Additionally, we recommend that you familiarize yourself with NetworKit through experimenting with the interactive IPython Notebook `NetworKit_UserGuide.ipynb` located in the folder `Doc/Notebooks`. The user guide also introduces a large portion of NetworKits functionality with usage examples. To display and work with these notebooks, you have to install jupyterhub and start a local notebook server from the terminal with:
::
	jupterhub --no-ssl

If you run into any problems with jupyterhub, head over to the `jupyterhub documentation <https://jupyterhub.readthedocs.io/en/latest/>`_ and make sure, you have the listed packages installed. If the notebook server starts as it is supposed to, your default browser should open a web interface or you have to open it manually. Then you can add `NetworKit_UserGuide.ipynb` from the above mentioned location or browse to the location through the web interface.

To show plots within the notebooks, place the following two lines at the beginning of your notebook:

.. code-block:: python

	%matplotlib
	matplotlib.pyplot as plt

**Note:** Instead of running jupyterhub, it may still be possible to run :code:`ipython3 notebook`. However, the notebook functionality of the ipython package is deprecated and has been moved to jupyterhub, which we strongly recommend.

Usage Example
~~~~~~~~~~~~~

Now that you are done installing NetworKit, you might want to try the following example:

	>>> from networkit import *
	>>> g = generators.HyperbolicGenerator(1e5).generate()
	>>> overview(g)
	Network Properties for:		    G#5
	nodes, edges			        100000, 302148
	directed?			            False
	weighted?			            False
	isolated nodes			        1859
	self-loops			            0
	density				            0.000060
	clustering coefficient		    0.718261
	min/max/avg degree		        0, 1045, 6.042960
	degree assortativity		    0.000725
	number of connected components	4237
	size of largest component	    77131 (77.13 %)
	>>> communities = community.detectCommunities(g, inspect=True)
	PLM(balanced,pc,turbo) detected communities in 0.3781468868255615 [s]
	solution properties:
	-------------------  -----------
	# communities        4468
	min community size      1
	max community size   1820
	avg. community size    22.3814
	modularity              0.989285
	-------------------  -----------
	>>>
