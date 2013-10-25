from NetworKit import *
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
					t = timer.elapsed

					mod = Modularity().getQuality(zeta, G)
					nc = zeta.numberOfClusters()


					row = [graphName, algoName, t, mod, nc]
					writer.writerow(row)
					print(row)

