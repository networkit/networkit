import csv
import os
try:
	import pandas
except ImportError:
	have_pandas = False
else:
	have_pandas = True
from networkit import *
from ..support import MissingDependencyError

def getFileList(directory):
	""" Get list of graph files in directory"""
	ls = []
	for (root, _, filenames) in os.walk(directory):
		for filename in filenames:
			ls.append(os.path.join(root, filename))
	return ls


def graphInfo(directory, fileEnding, outPath, fileformat=Format.METIS):
	if not have_pandas:
		raise MissingDependencyError("pandas")
	fileList = []
	data = []
	for (root, _, filenames) in os.walk(directory):
		for filename in filenames:
			if filename.endswith(fileEnding):
				fileList.append(os.path.join(root, filename))
	for file in fileList:
		print("reading ", file)
		G = readGraph(file, fileformat)
		print("analyzing ", file)
		graphName = os.path.split(file)[1].split(".")[0]
		data.append({"graph": graphName, "n": G.numberOfNodes(), "m": G.numberOfEdges(), "maxdeg": properties.degrees(G)[1], "comp": properties.numberOfComponents(G), "lcc": round(properties.clustering(G), 4)})
	print(data)
	info = pandas.DataFrame(data, columns=["graph", "n", "m", "maxdeg", "comp", "lcc"])
	info.to_csv(path_or_buf=outPath, sep="\t")
	return info


def communityDetectionBenchmark(graphPaths, outPath, algorithms, arglist=None, repeat=1, graphFormat=Format.METIS):
	"""
		Evaluate community detection algorithms on a collection of graphs and save benchmark data in .csv format
		:param	graphPaths	paths to graph files
		:param 	algorithms	list of algorithms
	"""
	if not arglist:
		arglist = [{} for i in range(len(algorithms))]

	# write results
	with open(outPath, 'w') as outFile:
		writer = csv.writer(outFile, delimiter='\t')
		header = ["graph", "algo", "time", "modularity"]
		writer.writerow(header)
		for graphPath in graphPaths:
			print("reading graph: {0}".format(graphPath))
			G = graphio.readGraph(graphPath, fileformat=graphFormat)
			graphName = os.path.basename(graphPath).split(".")[0]
			(n, m) = properties.size(G)
			for (algoClass, kwargs) in zip(algorithms, arglist):
				for i in range(repeat):
					algo = algoClass(G, **kwargs)
					algoName = algo.toString()
					print("evaluating {0} on {1}".format(algoName, graphName))
					timer = stopwatch.Timer()
					algo.run()
					zeta = algo.getPartition()
					timer.stop()
					time = timer.elapsed

					mod = community.Modularity().getQuality(zeta, G)
					# nc = zeta.numberOfSubsets()

					row = [graphName, algoName, time, mod]
					writer.writerow(row)
					print(row)


def testPLPThreshold(graphPaths, thresholdFactors, outPath, repeat=1):
	from NetworKit import community

	with open(outPath, 'w') as outFile:
		writer = csv.writer(outFile, delimiter='\t')
		header = ["graph", "algorithm","factor", "threshold", "time", "modularity"]
		writer.writerow(header)
		for graphPath in graphPaths:
			print("reading graph: {0}".format(graphPath))
			G = graphio.readGraph(graphPath)
			(n, m) = properties.size(G)
			graphName = os.path.basename(graphPath).split(".")[0]
			for factor in thresholdFactors:
				theta = int(n * factor)
				for i in range(repeat):
					timer = stopwatch.Timer()
					algo = community.PLP(updateThreshold=theta)
					algoName = algo.toString()
					zeta = algo.run(G)
					timer.stop()
					time = timer.elapsed

					mod = community.Modularity().getQuality(zeta, G)
					row = [graphName, "PLP", factor, theta,  time, mod]
					writer.writerow(row)
					print(row)


def testPLMDetailedScaling(graphPaths, threadSequence):
	data = {}
	for path in graphPaths:
		graphName = os.path.basename(path).split(".")[0]
		print("reading ", path)
		G = readGraph(path, Format.METIS)
		data[graphName] = {}
		for nThreads in threadSequence:
			setNumberOfThreads(nThreads)
			print("running on {0} at {1} threads".format(graphName, nThreads))
			plm = community.PLM(G, turbo=True, refine=True)
			plm.run()
			data[graphName][nThreads] = plm.getTiming()
	return data
