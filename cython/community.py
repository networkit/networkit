""" This module handles community detection, i.e. the discovery of densely connected groups in networks."""

from _NetworKit import Clustering, Coverage, Modularity, Clusterer, PLP, LPDegreeOrdered, PLM, MLPLM, ClusteringReader, ClusteringWriter, NodeStructuralRandMeasure, GraphStructuralRandMeasure


try:
	import tabulate
except ImportError:
	print(""" WARNING: module 'tabulate' not found, please install it to use the full functionality of NetworKit """)
import stopwatch

def detectCommunities(G, algo=None, inspect=False):
	""" Perform high-performance community detection on the graph.
		:param    G    the graph
		:param     algorithm    community detection algorithm instance
		:return communities (as type Clustering)
		"""
	if algo is None:
		algo = PLM()
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
	communitySizes = zeta.clusterSizes()
	mod = Modularity().getQuality(zeta, G)
	commProps = [
		["# communities", zeta.numberOfClusters()],
		["min community size", min(communitySizes)],
		["max community size", max(communitySizes)],
		["avg. community size", sum(communitySizes) / len(communitySizes)],
		["imbalance", zeta.getImbalance()],
		["modularity", mod],
	]
	print(tabulate.tabulate(commProps))
	

def evalCommunityDetection(algo, G):
	""" Evaluate a community detection algorithm """
	
	t = stopwatch.Timer()
	zeta = algo.run(G)
	t.stop()
	results = [
		["time [s]", t.elapsed],
		["# communities", zeta.numberOfClusters()],
		["modularity", Modularity().getQuality(zeta, G)]
	]
	print(tabulate.tabulate(results))


def readCommunities(path):
	""" Read a partition into communities from a file"""
	communities = ClusteringReader().read(path)
	print("read communities from: {0}".format(path))
	return communities


def writeCommunities(communities, path):
	""" Write a partition into communities to a file"""
	ClusteringWriter().write(communities, path)
	print("wrote communities to: {0}".format(path))


def compareCommunities(G, zeta1, zeta2):
	""" Compare the partitions with respect to several (dis)similarity measures"""
	pass # TODO