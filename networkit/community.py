""" This module handles community detection, i.e. the discovery of densely connected groups in networks."""

__author__ = "Christian Staudt"


from _NetworKit import Partition, Coverage, Modularity, CommunityDetector, PLP, LPDegreeOrdered, PLM, CNM, PartitionReader, PartitionWriter,\
	NodeStructuralRandMeasure, GraphStructuralRandMeasure, JaccardMeasure, NMIDistance,\
	EPP, EPPFactory, CommunityGraph, EdgeListPartitionReader, GraphClusteringTools, ClusteringGenerator, PartitionIntersection, HubDominance, CoreDecomposition, CutClustering

# local imports
#from .properties import CoreDecomposition, overview
from . import graph
from . import stopwatch

# external imports
import os
try:
	import tabulate
except ImportError:
	print(""" WARNING: module 'tabulate' not found, please install it to use the full functionality of NetworKit """)

def detectCommunities(G, algo=None, inspect=True):
	""" Perform high-performance community detection on the graph.
		:param    G    the graph
		:param     algorithm    community detection algorithm instance
		:return communities (as type Partition)
		"""
	if algo is None:
		algo = PLM(G, refine=False)
	t = stopwatch.Timer()
	algo.run()
	zeta = algo.getPartition()
	t.stop()
	print("{0} detected communities in {1} [s]".format(algo.toString(), t.elapsed))
	if inspect:
		print ("solution properties:")
		inspectCommunities(zeta, G)
	return zeta

def inspectCommunities(zeta, G):
	""" Display information about communities
		:param    zeta    communities
		:param    G        graph
	"""
	communitySizes = zeta.subsetSizes()
	mod = Modularity().getQuality(zeta, G)
	commProps = [
		["# communities", zeta.numberOfSubsets()],
		["min community size", min(communitySizes)],
		["max community size", max(communitySizes)],
		["avg. community size", sum(communitySizes) / len(communitySizes)],
		#["imbalance", zeta.getImbalance()],
		["modularity", mod],
	]
	print(tabulate.tabulate(commProps))


def communityGraph(G, zeta):
	""" Create a community graph, i.e. a graph in which one node represents a community and an edge represents the edges between communities, from a given graph and a community detection solution"""
	cg = CommunityGraph()
	cg.run(G, zeta)
	return cg.getGraph()


def evalCommunityDetection(algo, G):
	""" Evaluate a community detection algorithm """

	t = stopwatch.Timer()
	algo.run()
	zeta = algo.getPartition()
	t.stop()
	results = [
		["time [s]", t.elapsed],
		["# communities", zeta.numberOfSubsets()],
		["modularity", Modularity().getQuality(zeta, G)]
	]
	print(tabulate.tabulate(results))


def readCommunities(path, format="default"):
	""" Read a partition into communities from a file"""
	readers =  {"default": PartitionReader(),
		"edgelist-t1": EdgeListPartitionReader(1),
		"edgelist-t0": EdgeListPartitionReader(0),
		"edgelist-s1": EdgeListPartitionReader(1),
		"edgelist-s0": EdgeListPartitionReader(0),
		}
	# get reader
	try:
		reader = readers[format]#(**kwargs)
	except KeyError:
		raise Exception("unrecognized format: {0}".format(format))

	# get proper file path
	if ("~" in path):
		path = os.path.expanduser(path)
		print("path expanded to: {0}".format(path))
	# check if file path leads to a valid file
	if not os.path.isfile(path):
		raise IOError("{0} is not a file".format(path))
	else:
		with open(path, "r") as file:    # catch a wrong path before it crashes the interpreteri
			print("read communities from: {0}".format(path))
			communities = reader.read(path)
			return communities

	return None


def writeCommunities(communities, path):
	""" Write a partition into communities to a file"""
	PartitionWriter().write(communities, path)
	print("wrote communities to: {0}".format(path))


def compareCommunities(G, zeta1, zeta2):
	""" Compare the partitions with respect to several (dis)similarity measures"""
	raise NotImplementedError("TODO:")


def kCoreCommunityDetection(G, k, algo=None, inspect=True):
	""" Perform community detection on the k-core of the graph, which possibly
		reduces computation time and enhances the result.
		:param    G    the graph
		:param		k 	k as in k-core
		:param     algorithm    community detection algorithm instance
		:return communities (as type Partition)
		"""
	coreDec = CoreDecomposition(G)
	coreDec.run()

	cores = coreDec.cores()
	try:
		kCore = cores[k]
	except IndexError:
		raise Error("There is no core for the specified k")

	C = graph.Subgraph().fromNodes(G, kCore)	# FIXME: node indices are not preserved

	#properties.overview(C)

	return detectCommunities(C, algo, inspect)
