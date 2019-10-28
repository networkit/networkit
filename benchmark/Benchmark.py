"""
The NetworKit benchmark

This module implements a comprehensive benchmark of NetworKit's analytics algorithms



"""


import sys
import warnings
import math
import os
import numpy

from time import gmtime, strftime
import signal
import networkit

from util import *
from networkit.support import MissingDependencyError
import nk

try:
	import pandas
except ImportError:
	have_pandas = False
else:
	have_pandas = True
try:
	import matplotlib.pyplot as plt
except ImportError:
	have_plt = False
else:
	have_plt = True
try:
	import seaborn
except ImportError:
	have_seaborn = False
else:
	have_seaborn = True

try:
	import networkx
	import nx
except ImportError:
	have_nx = False
else:
	have_nx = True

try:
	import igraph
	import ig
except ImportError:
	have_ig = False
else:
	have_ig = True

try:
	import graph_tool
	import gt
except ImportError:
	have_gt = False
else:
	have_gt = True

# helper function


def averageRuns(df, groupby=["graph"]):
	""" Average running time, modularity, edges per second and number of clusters over multiple runs"""
	df = df.groupby(groupby, as_index=False).mean()
	df = df.sort("m", ascending=True)	# sort by graph size
	return df


def graphMeta(graphNames, graphDir, fileEnding=".gml.graph", graphFormat=networkit.Format.GML):
	if not have_pandas:
		raise MissingDependencyError("pandas")
	meta = []
	for name in graphNames:
		info("loading {name}".format(**locals()))
		G = networkit.readGraph(os.path.join(graphDir, "{0}{1}".format(name, fileEnding)), graphFormat)
		(n, m) = G.size()
		meta.append({"name" : name, "n" : n, "m" : m})
	info("done")
	return pandas.DataFrame(meta, columns=["name", "n", "m"])





def saveData(df, name):
	df.to_csv(os.path.join(bench.dataPath, "{0}.csv".format(name)), sep="\t")

# timeout

class Timeout(Exception):
	pass

def timeoutHandler(signum, frame):
	error("timeout")
	raise Timeout()




class bFail:
	name = "Fail"

	def run(self, G):
		raise Exception("FAIL!")




# Plots

## plot settings
if not have_seaborn:
	raise MissingDependencyError("seaborn")
seaborn.set_style("whitegrid")

### Colors
lightred = seaborn.xkcd_rgb["red"]
darkred = seaborn.xkcd_rgb["crimson"]
green = seaborn.xkcd_rgb["teal"]
orange = seaborn.xkcd_rgb["bright orange"]

# plot functions

def timePlot(data, size=(6,3)):
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	pos = numpy.arange(len(data))+.5	# the bar centers on the y axis
	labels = list(data["graph"])
	plt.figure(figsize=size)
	plt.xscale("symlog")
	plt.barh(pos, data["time"], align='center', height=0.25, color=lightred)	# notice the 'height' argument
	plt.yticks(pos, labels)
	plt.gca().xaxis.set_minor_locator(plt.LogLocator(subs=[0,1,2,3,4,5,6,7,8,9,10]))
	#gca().xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
	plt.xlabel("time [s]")
	plt.grid(True)


def epsPlot(data, size=(6,3)):
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	pos = numpy.arange(len(data))+.5	# the bar centers on the y axis
	labels = list(data["graph"])
	plt.figure(figsize=size)
	plt.xscale("log")
	plt.barh(pos, data["edges/s"], align='center', height=0.25, color=green)	# notice the 'height' argument
	plt.yticks(pos, labels)
	plt.gca().xaxis.set_minor_locator(plt.LogLocator(subs=[0,1,2,3,4,5,6,7,8,9,10]))
	#gca().xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
	plt.xlabel("edges / s")
	plt.grid(True)



def graphPlot(data, size=None):
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	pos = arange(len(data))	# the bar centers on the y axis
	plt.figure(figsize=(5,3.3))
	sizes = data["m"]
	labels = list(data["graph"])
	#barh(pos, data["n"], log=True, align='center', height=0.4, color="darkgrey")	# notice the 'height' argument
	plt.barh(pos, sizes, log=True, align='center', height=0.4, color="lightblue", label="m")	# notice the 'height' argument
	plt.yticks(pos, labels)
	plt.xlabel("number of edges")
	plt.grid(True)


def barPlot(data, labels, x_label="", y_label="", size=None, color="b", transparency=1.0, scale="linear"):
	if not have_plt:
		raise MissingDependencyError("matplotlib")
	pos = numpy.arange(len(data))+.5	# the bar centers on the y axis
	plt.figure(figsize=size)
	plt.xscale(scale)
	plt.barh(pos, data, align='center', height=0.25, color=color, alpha=transparency)	# notice the 'height' argumen
	plt.yticks(pos, labels)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.grid(True)


class Bench:

	"""
	A Bench objects represents a benchmarking session. It is responsible
	for running the benchmarks, collecting and storing the resulting data.

	"""

	def __init__(self, graphDir, defaultGraphs, outDir, save=True, nRuns=5, cacheGraphs=False, timeout=None, fileEnding="gml.graph", graphFormat=networkit.Format.GML, indexEdges=False):
		self.graphFormat = graphFormat
		self.indexEdges = indexEdges
		self.fileEnding = fileEnding
		self.defaultGraphs = defaultGraphs
		self.nRuns = nRuns  # default number of runs for each algo
		self.graphDir = graphDir
		# store result data of benchmarks
		self.data = []  # list of dictionaries (rows) with keys representing table columns
		self.loadTimes = []
		self.outDir = outDir
		self.save = save  # store data frames on disk if true
		# create output directory if it does not exist
		if self.save:
			if not os.path.isdir(self.outDir):
				os.mkdir(self.outDir)
			self.outDataDir = os.path.join(self.outDir, "data")
			self.plotDir = os.path.join(self.outDir, "plots")
			if not os.path.isdir(self.outDataDir):
				os.mkdir(self.outDataDir)
			if not os.path.isdir(self.plotDir):
				os.mkdir(self.plotDir)
			# log file
			self.logPath = os.path.join(self.outDir, "log.txt")
		# graph cache
		self.cacheGraphs = cacheGraphs
		self.graphCache = {}
		#timeout
		self.timeout = timeout

	def getGraph(self, graphName, algo):
		"""" Get the graph from disk or from in-memory cache"""
		key = algo.framework + graphName
		if key in self.graphCache:
			self.info("getting {0} from cache".format(graphName))
			G = self.graphCache[key]
			return G
		else:
			graphPath = os.path.join(self.graphDir, "{0}.{1}".format(graphName, self.fileEnding))
			self.info("loading {0} from {1}".format(graphName, graphPath))
			with Timer() as t:
				G = algo.loadGraph(graphPath, self.graphFormat)
				if self.indexEdges:
					G.indexEdges()
			self.info("\t took {1} s".format(graphName, t.elapsed))
			eps = self.getSize(G)[1] / t.elapsed
			self.loadTimes.append({"framework" : algo.framework, "graph" : graphName, "time": t.elapsed, "edges/s": eps})
			if self.cacheGraphs:
				self.graphCache[key] = G
			return G

	def getSize(self, G):
		if not have_gt:
			raise MissingDependencyError("graph_tool")
		if not have_ig:
			raise MissingDependencyError("igraph")
		if isinstance(G, networkit.Graph):
			return G.size()
		elif isinstance(G, igraph.Graph):
			return (G.vcount(), G.ecount())
		elif isinstance(G, graph_tool.Graph):
			return (G.num_vertices(), G.num_edges())
		else:
			raise NotImplementedError("cannot get graph size - unknown object type")

	def clearCache(self):
		""" Delete all stored graphs to free memory """
		del self.graphCache
		self.graphCache = {}

	def algoBenchmark(self, algo, graphs=None, nRuns=None, timeout=None):
		""" Run a kernel represented by an algorithm benchmark object """
		# set the defaults
		if nRuns is None:
			nRuns = self.nRuns  # lets argument override the default nRuns
		if graphs is None:
			graphs = self.defaultGraphs
		if timeout is None:
			timeout = self.timeout

		self.info("benchmarking {algo.name}".format(**locals()))
		#table = []  # list of dictionaries, to be converted to a DataFrame

		for graphName in graphs:
			try:
				G = self.getGraph(graphName, algo)
				(n, m) = self.getSize(G)
				try:
					self.info("running {algo.name} ({algo.framework}) {nRuns}x on {graphName}".format(**locals()))
					for i in range(nRuns):
						row = {}	# benchmark data row
						row["algorithm"] = algo.name
						row["framework"] = algo.framework
						row["graph"] = graphName
						row["m"] = m
						try: # timeout
							if timeout:
								signal.signal(signal.SIGALRM, timeoutHandler)
								signal.alarm(int(timeout * 60))  # timeout in seconds
							with Timer() as t:
								result = algo.run(G)
							self.info("\t took {0} s".format(t.elapsed))
							# store data
							row["time"] = t.elapsed
							row["edges/s"] =  m / t.elapsed  # calculate edges per second
							row["result"] = result
						except Timeout as tx:
							self.error("{algo.name} timed out after {timeout} minutes".format(**locals()))
							row["time"] = None
							row["edges/s"] = None
							row["result"] = None
						finally:
							self.data.append(row)
							signal.alarm(int(1e9))	# in any case, cancel the timeout alarm by setting it to a ridiculously high time
				except Exception as ex:
					self.error("algorithm {algo.name} failed with exception: {ex}".format(**locals()))
			except Exception as ex:
				self.error("loading graph {graphName} failed with exception: {ex}".format(**locals()))

		#df = pandas.DataFrame(table)
		#df.sort("m")	# sort by number of edges
		#self.data[algo.name] = df
		# store data frame on disk
		#if self.save:
		#	df.to_csv(os.path.join(self.outDataDir, "{algo.name}.csv".format(**locals())))
		self.log("Done")


	def log(self, message):
		""" Write a message to the logfile"""
		with open(self.logPath, "a") as logfile:
			logfile.write("{0}: {1}\n".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()), message))

	def info(self, message):
		print(message)
		sys.stdout.flush()
		self.log(message)

	def error(self, message):
		print(message)
		sys.stdout.flush()
		self.log(message)

	def debug(self, message):
		pass
		# print(message)


	def timePlot(self, algoName):
		if not have_plt:
			raise MissingDependencyError("matplotlib")
		timePlot(averageRuns(self.data[algoName]))
		if self.save:
			plt.savefig(os.path.join(self.plotDir, "{algoName}-time.pdf".format(**locals())), bbox_inches="tight")

	def epsPlot(self, algoName):
		if not have_plt:
			raise MissingDependencyError("matplotlib")
		epsPlot(averageRuns(self.data[algoName]))
		if self.save:
			plt.savefig(os.path.join(self.plotDir, "{algoName}-eps.pdf".format(**locals())), bbox_inches="tight")

	def graphPlot(self):
		if not have_plt:
			raise MissingDependencyError("matplotlib")
		epsPlot(averageRuns(self.data))
		if self.save:
			plt.savefig(os.path.join(self.plotDir, "graphs.pdf".format(**locals())), bbox_inches="tight")

	def plot(self, algoName):
		print(algoName)
		self.timePlot(algoName)
		self.epsPlot(algoName)

	def finalize(self):
		if not have_pandas:
			raise MissingDependencyError("pandas")
		self.dataFrame = pandas.DataFrame(self.data)
		if self.save:
			self.dataFrame.to_csv(os.path.join(self.outDataDir, "data.csv".format(**locals())))

	def plotSummary2(self, figsize=None, groupby="framework", palette="Greens_d"):
		""" Plot a summary of algorithm performances"""
		if not have_plt:
			raise MissingDependencyError("matplotlib")
		if not have_seaborn:
			raise MissingDependencyError("seaborn")
		if figsize:
			plt.figure(figsize=figsize)
		plt.gca().xaxis.get_major_formatter().set_powerlimits((3, 3))
		plt.xscale("log")
		plt.xlabel("edges/s")
		ax = seaborn.boxplot(y="algorithm", x="edges/s", hue=groupby, data=self.dataFrame, linewidth=1, width=.5, palette=palette)
		if self.save:
			plt.savefig(os.path.join(self.plotDir, "epsSummary.pdf".format(**locals())), bbox_inches="tight")

	def plotSummary(self, algoNames=None, figsize=None):
		""" Plot a summary of algorithm performances"""
		if not have_plt:
			raise MissingDependencyError("matplotlib")
		if not have_pandas:
			raise MissingDependencyError("pandas")
		if not have_seaborn:
			raise MissingDependencyError("seaborn")
		if algoNames is None:
			algoNames = list(self.data.keys())
		epsSummary = pandas.DataFrame()
		for (algoName, algoData) in self.data.items():
			if algoName in algoNames:
				epsSummary[algoName] = pandas.Series(self.data[algoName]["edges/s"])

		# data frame
		self.epsSummary = epsSummary
		self.epsSummary = self.epsSummary.reindex_axis(sorted(self.epsSummary.columns), axis=1)
		if self.save:
			self.epsSummary.to_csv(os.path.join(self.outDataDir, "epsSummary.csv".format(**locals())))
		# plot
		if figsize:
			plt.figure(figsize=figsize)
		plt.gca().xaxis.get_major_formatter().set_powerlimits((3, 3))
		plt.xscale("log")
		plt.xlabel("edges/s")
		seaborn.boxplot(self.epsSummary, order=self.epsSummary.columns, linewidth=1.5, width=.25, color=green, vert=False)
		if self.save:
			plt.savefig(os.path.join(self.plotDir, "epsSummary.pdf".format(**locals())), bbox_inches="tight")

	def plotAll(self):
		for key in self.data.keys():
			self.plot(key)


	def getLoadTimes(self):
		""" Get input times for graphs"""
		if not have_pandas:
			raise MissingDependencyError("pandas")
		return pandas.DataFrame(self.loadTimes)



	def generatorBenchmark(self, generator, argtuples, nRuns=None, timeout=None):
		""" Run a kernel represented by an algorithm benchmark object """
		if not have_pandas:
			raise MissingDependencyError("pandas")
		# set the defaults
		if nRuns is None:
			nRuns = self.nRuns  # lets argument override the default nRuns
		if timeout is None:
			timeout = self.timeout

		genName = str(generator)

		self.info("benchmarking {genName}".format(**locals()))
		table = []  # list of dictionaries, to be converted to a DataFrame

		for param in argtuples:
			try:
				gen = generator(*param)
				row = {}	# benchmark data row
				row["graph"] = str(param)
				try:
					self.info("running {genName} {nRuns} times".format(**locals()))
					for i in range(nRuns):
						row["algorithm"] = genName
						try: # timeout
							result = None
							if timeout:
								signal.signal(signal.SIGALRM, timeoutHandler)
								signal.alarm(int(timeout * 60))  # timeout in seconds
							with Timer() as t:
								result = gen.generate()
							self.debug("took {0} s".format(t.elapsed))
							# store data
							row["m"] = result.numberOfNodes()
							row["time"] = t.elapsed
							row["result"] = result
						except Timeout as tx:
							self.error("{genName} timed out after {timeout} minutes".format(**locals()))
							row["time"] = None
							row["result"] = None
						finally:
							table.append(row)
							signal.alarm(int(1e9))	# in any case, cancel the timeout alarm by setting it to a ridiculously high time
				except Exception as ex:
					self.error("generator {genName} failed with exception: {ex}".format(**locals()))
			except Exception as ex:
				self.error("initializing generator {genName} failed with exception: {ex}".format(**locals()))
		if not have_pandas:
			raise MissingDependencyError("pandas")
		df = pandas.DataFrame(table)
		self.data[genName] = df
		# store data frame on disk
		if self.save:
			df.to_csv(os.path.join(self.outDataDir, "{genName}.csv".format(**locals())))



# - generators
# 	- Erdös-Renyi (generators.ErdosRenyiGenerator)
# 	- Barabasi-Albert (generators.BarabasiAlbertGenerator)
# 	- Chung-Lu (generators.ChungLuGenerator)
# 	- RMAT (generators.RmatGenerator)
