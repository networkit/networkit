""" This module handles community detection, i.e. the discovery of densely connected groups in networks."""

from _NetworKit import Partition, Coverage, Modularity, CommunityDetector, PLP, LPDegreeOrdered, PLM, CNM, PartitionReader, PartitionWriter, NodeStructuralRandMeasure, GraphStructuralRandMeasure, JaccardMeasure, EPP, EPPFactory, CommunityGraph, EdgeListPartitionReader
import os

try:
	import tabulate
except ImportError:
	print(""" WARNING: module 'tabulate' not found, please install it to use the full functionality of NetworKit """)
import stopwatch

def detectCommunities(G, algo=None, inspect=True):
	""" Perform high-performance community detection on the graph.
		:param    G    the graph
		:param     algorithm    community detection algorithm instance
		:return communities (as type Clustering)
		"""
	if algo is None:
		algo = PLM(refine=False)
	t = stopwatch.Timer()
	zeta = algo.run(G)
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
	zeta = algo.run(G)
	t.stop()
	results = [
		["time [s]", t.elapsed],
		["# communities", zeta.numberOfSubsets()],
		["modularity", Modularity().getQuality(zeta, G)]
	]
	print(tabulate.tabulate(results))


def readCommunities(path, format="partition"):
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
	pass # TODO
