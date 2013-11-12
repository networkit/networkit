""" 
NetworKit - an interactive toolkit for high-performance network analysis
"""


__author__ = "Christian L. Staudt (christian.staudt @ kit.edu)"
__copyright__ = "Copyright (c) 2013 Christian Staudt"
__license__ = "MIT License"
__version__ = "2.0" 


# standard library modules


# non standard library modules
try:
	import networkx as nx
except ImportError:
	print("""WARNING: module 'networkx' not installed, which is required by some
						functions.""")

# faulthandler prints debug info on errors like segfaults
try:
	import faulthandler
	faulthandler.enable()
except ImportError:
	print("""" WARNING: module 'faulthandler' not found. If installed, it will 
		provide additional debug info on errors like segmentation faults.""")

# local modules
import stopwatch


# import from extension module into NetworKit namespace
from _NetworKit import configureLogging, currentLogLevel

# NetworKit submodules
import graph
from graph import Graph 
import graphio
import community
import generators
import properties


#-------- Setup ---------- #

def setup():
	""" This function is run once on module import to configure initial settings """
	configureLogging("ERROR")    # set default loglevel to error
	

setup() # here the setup function is called once on import


#--------- NetworKit Python Shell functions ----------------#

def readGraph(path, format=None):
	"""    Read graph and return a NetworKit::Graph"""
	# TODO: detect file format by looking at the file content
	if format is None:    # look at file extension
		if path.endswith(".graph") or path.endswith(".metis") or path.endswith(".dimacs"):
			reader = graphio.METISGraphReader()
		else:
			raise Exception("unknown graph file format")
	else:
		if (format is "metis"):
			reader = graphio.METISGraphReader()
		elif (format is "edgelist"):
			reader = graphio.EdgeListIO('\t', 1)

	with open(path, "r") as file:    # catch a wrong path before it crashes the interpreter
		G = reader.read(path)
		return G

	return None

########  CONVERSION ########

def nx2nk(nxG, weightAttr=None):
	""" 
	Convert a networkx.Graph to a NetworKit.Graph
		:param weightAttr: the edge attribute which should be treated as the edge weight
	 """
	# TODO: consider weights
	n = nxG.number_of_nodes()
	nkG = Graph(n)
	
	if weightAttr is not None:
		nkG.markAsWeighted()
		for (u, v) in nxG.edges():
			w = nxG.edge[u][v][weightAttr]
			nkG.addEdge(u, v, w)
	else:
		for (u, v) in nxG.edges():
			nkG.addEdge(u, v)
	
	return nkG


def nk2nx(nkG):
	""" Convert a NetworKit.Graph to a networkx.Graph """
	nxG = nx.Graph()
	if nkG.isMarkedAsWeighted():
		for (u, v) in nkG.edges():
			nxG.add_edge(u, v, weight=nkG.weight(u, v))
	else:
		for (u, v) in nkG.edges():
			nxG.add_edge(u, v)
	return nxG
	




# working with attributes

def retrieveAttributes(nodes, attributes):
	"""
	For a given collection of nodes and a map (or vector) indexed by node ids
	return a map from node to attribute.
	"""
	attrs = dict()
	for u in nodes:
		attrs[u] = attributes[u]
	return attrs
	
	
	
