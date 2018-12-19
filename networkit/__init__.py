"""
NetworKit -- an interactive tool suite for high-performance network analysis.

NetworKit is an open-source software package for high-performance
analysis of large complex networks. Complex networks are equally
attractive and challenging targets for data mining, and novel
algorithmic solutions, including parallelization, are required to handle
data sets containing billions of connections. Our goal for NetworKit is
to package results of our algorithm engineering efforts and put them
into the hands of domain experts. NetworKit is a hybrid combining the
performance of kernels written in C++ with a convenient Python
frontend. The package targets shared-memory platforms with OpenMP support.
The current feature set includes various analytics kernels such as connected components, diameter, clustering
coefficients, community detection, k-core decomposition, degree assortativity
and multiple centrality indices, as well as a
collection of graph generators. Scaling to massive networks is enabled by techniques such as parallel and sampling-based approximation algorithms.
NetworKit is geared towards large networks and satisfies three important
criteria: High performance, interactive workflows
and integration into the Python ecosystem of tools for data analysis and
scientific computation.


Usage examples can be found on http://nbviewer.ipython.org/urls/networkit.github.io/data/uploads/docs/NetworKit_UserGuide.ipynb
"""

__author__ = "Christian Staudt"
__copyright__ = "Copyright (c) 2014 Christan Staudt"
__credits__ = ["Eugenio Angriman", "Lukas Barth", "Miriam Beddig", "Elisabetta Bergamini", "Stefan Bertsch", "Pratistha Bhattarai", "Andreas Bilke", "Simon Bischof", \
	"Guido Brückner", "Mark Erb",  "Kolja Esders", "Patrick Flick", "Michael Hamann", "Lukas Hartmann", "Daniel Hoske", "Gerd Lindner", "Moritz v. Looz", "Yassine Marrakchi", "Henning Meyerhenke", \
	"Manuel Penschuck", "Marcel Radermacher", "Klara Reichard", "Marvin Ritter", "Aleksejs Sazonovs", "Hung Tran", "Florian Weber", "Michael Wegner", "Jörg Weisbarth", "Kolja Esders"]
__license__ = "MIT"
__version__ = "5.0"

# standard library modules
import csv
import os
import logging
import sys

# pyplot might not be available in non-interactive environments.
try:
	import matplotlib.pyplot as _pyplot
except ImportError as importError:
	_pyplot = None

# local imports
from . import stopwatch
from . import graph
from . import graphio
from . import community
from . import centrality
from . import generators
from . import structures
from . import engineering
from . import distance
from . import components
from . import dynamic
from . import gephi
from . import partitioning
from . import coloring
from . import workflows
from . import flow
from . import sparsification
from . import scd
from . import clique
from . import globals
from . import linkprediction
from . import correlation
from . import matching
from . import coarsening
from . import simulation
from . import stats
from . import sampling
from . import viz
from . import randomization

if _pyplot is not None:
	from . import plot
	from .profiling import profiling
else:
	print("WARNING: pyplot is not available; some functionality is disabled", file=sys.stderr)

try:
	from . import viztasks
except ImportError as importError:
	print("WARNING: some dependencies are not satisfied"
			" which are needed to use the 'viztasks' submodule", file=sys.stderr)

#--------- Top Level Classes and Functions ----------------#
#

# Some functions and classes should be directly available from the top module

# TODO: introduce settings module

# extension imports
from _NetworKit import getLogLevel, setLogLevel, setPrintLocation, enableNestedParallelism, setNumberOfThreads, getCurrentNumberOfThreads, getMaxNumberOfThreads, none, setSeed

# local imports into the top namespace
from .graph import Graph
from .structures import Partition, Cover
from .graphio import readGraph, writeGraph, readGraphs, Format


def overview(G):
	"""
		This function collects some basic information about the given graph and prints it to the terminal.
	"""
	n = G.numberOfNodes()
	degrees = centrality.DegreeCentrality(
		G, ignoreSelfLoops=G.numberOfSelfLoops() == 0).run().scores()
	numSelfLoops = G.numberOfSelfLoops()

	def getIsolatedNodes(degrees):
		sequence = sorted(degrees)
		i = 0
		nIsolated = 0
		while i < len(sequence) and sequence[i] == 0:
			nIsolated += 1
			i += 1
		return nIsolated

	def getClusteringCoefficient(G):
		lcc = centrality.LocalClusteringCoefficient(G, True).run().scores()
		return sum(lcc) / n

	def getComponentPartition(G):
		if G.isDirected():
			cc = components.StronglyConnectedComponents(G).run()
		else:
			cc = components.ConnectedComponents(G).run()
		return cc.getPartition()

	print("Network Properties for:\t\t{}".format(G.getName()))
	print("nodes, edges\t\t\t{}, {}".format(n, G.numberOfEdges()))
	print("directed?\t\t\t{}".format("True" if G.isDirected() else "False"))
	print("weighted?\t\t\t{}".format("True" if G.isWeighted() else "False"))
	print("isolated nodes\t\t\t{}".format(getIsolatedNodes(degrees)))
	print("self-loops\t\t\t{}".format(numSelfLoops))
	print("density\t\t\t\t{:.6f}".format(G.density()))
	if numSelfLoops == 0 and not G.isDirected():
		print("clustering coefficient\t\t{:.6f}".format(
			getClusteringCoefficient(G)))
	print("min/max/avg degree\t\t{:d}, {:d}, {:.6f}".format(
		int(min(degrees)), int(max(degrees)),
		sum(degrees) / n))
	print("degree assortativity\t\t{:.6f}".format(
		correlation.Assortativity(G, degrees).run().getCoefficient()))
	cp = getComponentPartition(G)
	lcs = max(cp.subsetSizes())
	print("number of connected components\t{}".format(cp.numberOfSubsets()))
	print("size of largest component\t{} ({:.2f} %)".format(
		lcs, 100 * lcs / n))


#-------- Setup ---------- #


def setup():
	""" This function is run once on module import to configure initial settings """
	setLogLevel("ERROR")  # set default loglevel for C++ code
	setPrintLocation(True)
	enableNestedParallelism()  # enable nested parallelism
	logging.basicConfig(
		level=logging.INFO)  # set default loglevel for Python code


setup()  # here the setup function is called once on import

# in general, no implementations here
