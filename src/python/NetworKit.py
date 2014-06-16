"""
NetworKit - an interactive toolkit for high-performance network analysis
"""


__author__ = "Christian L. Staudt (christian.staudt @ kit.edu)"
__copyright__ = "Copyright (c) 2013 Christian Staudt"
__license__ = "MIT License"
__version__ = "3.2"


# standard library modules
import csv
import os
import logging

# local modules
import stopwatch

# NetworKit submodules
import graph
import graphio
import community
import centrality
import generators
import properties
import engineering
import toolbox
import dynamic
import gephi

try:
	import viztools
	import viztasks
	import algebraic
except ImportError as importError:
	print("""WARNING: some dependencies are not satisfied which are needed to use the
		'viztools' submodule""")
	print(importError)




#--------- Top Level Classes and Functions ----------------#
#

# Some functions and classes should be directly available from the top module

# TODO: introduce settings module

from _NetworKit import getLogLevel, setLogLevel, setPrintLocation, enableNestedParallelism, setNumberOfThreads, getCurrentNumberOfThreads, getMaxNumberOfThreads
from graph import Graph
from graphio import readGraph, writeGraph, formats
from nxadapter import nk2nx, nx2nk
from toolbox import batch
from community import detectCommunities


#-------- Setup ---------- #

def setup():
	""" This function is run once on module import to configure initial settings """
	setLogLevel("ERROR")    # set default loglevel for C++ code
	setPrintLocation(True)
	enableNestedParallelism()	# enable nested parallelism
	logging.basicConfig(level=logging.INFO)	# set default loglevel for Python code



setup() # here the setup function is called once on import


# in general, no implementations here
