"""
an interactive toolkit for high-performance network analysis
usage examples can be found on http://nbviewer.ipython.org/urls/networkit.iti.kit.edu/data/uploads/docs/NetworKit_UserGuide.ipynb
"""

__author__ = "Christian Staudt"
__copyright__ = "Copyright (c) 2014 Christan Staudt"
__license__ = "MIT"
__version__ = "3.4.1"


# standard library modules
import csv
import os
import logging
import sys

# local imports
from . import stopwatch
from . import graph
from . import graphio
from . import community
from . import centrality
from . import generators
from . import properties
from . import structures
from . import engineering
from . import dynamic
from . import gephi
from . import partitioning
from . import coloring
from . import workflows
from . import flow
from . import distmeasures
from . import plot
from . import scd

try:
	from . import viztools
	from . import viztasks
	from . import algebraic
except ImportError as importError:
	print("""WARNING: some dependencies are not satisfied which are needed to use the
		'viztools' submodule""")
	print(importError)




#--------- Top Level Classes and Functions ----------------#
#

# Some functions and classes should be directly available from the top module

# TODO: introduce settings module

# extension imports
from _NetworKit import getLogLevel, setLogLevel, setPrintLocation, enableNestedParallelism, setNumberOfThreads, getCurrentNumberOfThreads, getMaxNumberOfThreads, none

# local imports
from .graph import Graph
#try:
from .graphio import readGraph, writeGraph, Format
#except ImportError:
#	from _graphio33 import readGraph, writeGraph, Format
from .nxadapter import nk2nx, nx2nk
from .workflows import batch
from .community import detectCommunities


#-------- Setup ---------- #

def setup():
	""" This function is run once on module import to configure initial settings """
	setLogLevel("ERROR")    # set default loglevel for C++ code
	setPrintLocation(True)
	enableNestedParallelism()	# enable nested parallelism
	logging.basicConfig(level=logging.INFO)	# set default loglevel for Python code



setup() # here the setup function is called once on import


# in general, no implementations here
