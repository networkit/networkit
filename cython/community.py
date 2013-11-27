""" This module handles community detection, i.e. the discovery of densely connected groups in networks."""

from _NetworKit import Clustering, Coverage, Modularity, Clusterer, PLP, LPDegreeOrdered, PLM, MLPLM, \
ClusteringReader, ClusteringWriter, NodeStructuralRandMeasure, GraphStructuralRandMeasure, EPP, EPPFactory

try:
	import tabulate
except ImportError:
	print(""" WARNING: module 'tabulate' not found, please install it to use the full functionality of NetworKit """)
import stopwatch

def detectCommunities(G, algorithm=None, inspect=True):
	""" Perform high-performance community detection on the graph.
		:param    G    the graph
		:param     algorithm    community detection algorithm instance
		:return communities (as type Clustering)
		"""
	if algorithm is None:
		algorithm = PLM()
	t = stopwatch.Timer()
	zeta = algorithm.run(G)
	t.stop()
	print("{0} detected communities in {1} [s]".format(algorithm.toString(), t.elapsed))
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


import stopwatch
import csv
import os


def getFileList(directory):
	""" Get list of graph files in directory"""
	ls = []
	for (root, _, filenames) in os.walk(directory):
		for filename in filenames:
			ls.append(os.path.join(root, filename))
	return ls

def communityDetectionBenchmark(graphPaths, algorithms, outPath, repeat=1):
	"""
		Evaluate community detection algorithms on a collection of graphs and save benchmark data in .csv format
		:param	graphPaths	paths to graph files
		:param 	algorithms	list of algorithms
	"""

	# write results
	with open(outPath, 'w') as outFile:
		writer = csv.writer(outFile, delimiter=';')
		for graphPath in graphPaths:
			print("reading graph: {0}".format(graphPath))
			G = readGraph(graphPath)
			graphName = os.path.basename(graphPath).split(".")[0]
			(n, m) = nm(G)
			for algo in algorithms:
				algoName = algo.toString()
				for i in range(repeat):
					print("evaluating {0} on {1}".format(algoName, graphName))
					timer = stopwatch.Timer()
					zeta = algo.run(G)
					timer.stop()
					time = timer.elapsed

					mod = Modularity().getQuality(zeta, G)
					#nc = zeta.numberOfClusters()


					row = [graphName, algoName, time, mod]
					writer.writerow(row)
					print(row)


