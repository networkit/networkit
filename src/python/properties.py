# NetworKit native classes and functions
from _NetworKit import GraphProperties, ConnectedComponents, ClusteringCoefficient, Diameter, Eccentricity, CoreDecomposition

# other submodules
import community
import termgraph
import auxiliary
import nxadapter
import powerlaw

# other modules
import textwrap
import collections
import math
import logging

try:
	import networkx as nx
except ImportError:
	logging.warning("""WARNING: module 'networkx' not installed, which is required by some
						functions.""")
try:
	import tabulate
except ImportError:
	logging.warning("""WARNING: module 'tabulate' not installed, which is required by some
						functions. In order to install 'tabulate', Python 3.3 is required""")

try:
	from scipy import stats
except ImportError:
	logging.warning("""WARNING: module 'scipy' not installed, which is required by some
						functions.""")


########  PROPERTIES ########

def size(G):
	""" Return number of nodes and number of edges"""
	n = G.numberOfNodes()
	m = G.numberOfEdges()
	return (n, m)

def degrees(G):
	""" Return min/max/avg degree"""
	minMaxDeg = GraphProperties.minMaxDegree(G)
	avgDeg = GraphProperties.averageDegree(G)
	return (minMaxDeg[0], minMaxDeg[1], avgDeg)

def degreeDistribution(G):
	""" Return the degree distribution of a graph"""
	return GraphProperties.degreeDistribution(G)


def density(G):
	""" Return the density of the graph"""
	(n, m) = size(G)
	return (2 * m) / (n * (n-1))

def components(G):
	""" Find and analyze detected components.
		Returns the number of components and the sizes
		of each component. For more details use the
		ConnectedComponents class.
	"""
	logging.info("[...] finding connected components....")
	cc = ConnectedComponents(G)
	cc.run()
	components = cc.getPartition()
	nComponents = components.numberOfSubsets()
	componentSizes = components.subsetSizeMap()
	return (nComponents, componentSizes)

def numberOfComponents(G):
	""" Find and number of components """
	logging.info("[...] finding connected components....")
	cc = ConnectedComponents()
	cc.run(G)
	nComponents = cc.numberOfComponents()
	return nComponents

def clustering(G, error=0.01):
	"""
		Returns approximate average local clustering coefficient
		The maximum error can be given as a parameter and determines
		the number of samples taken.

		for details see:
			Schank, Wagner: Approximating Clustering Coefficient and Transitivity
	"""
	nSamples = math.ceil(math.log(10) / (error**2)) # fixed confidence of 90%
	logging.info("taking {0} samples".format(nSamples))
	return ClusteringCoefficient().approxAvgLocal(G, nSamples)


def powerLawExponent(G, dd=None):
	""" Supposing that the degreee distribution fits a power law,
		get the exponent (gamma).
	"""
	if not dd:
		dd = degreeDistribution(G)
	fit = powerlaw.Fit(dd)
	return fit.alpha

def powerLawFit(G, dd=None):
	""" Check if a power law is a good fit for the degree distribution.
	"""
	if not dd:
		dd = degreeDistribution(G)
	fit = powerlaw.Fit(dd)
	R, p = fit.distribution_compare("power_law", "exponential", normalized_ratio=True)
	return ((R > 0), R)


# def powerLawExponent_(G):
#     """ Estimate power law exponent as a linear regression line through a log-log plot
#     of the degree distribution"""
#     def logarithmic(data):
#     	return [(math.log(x) if x != 0 else 0) for x in data]

#     dd = degreeDistribution(G)
#     degrange = range(len(dd))
#     (slope,_,_,_,_) = stats.linregress(logarithmic(degrange), logarithmic(dd))
#     gamma = -1 * slope
#     return gamma


def degreeAssortativity(G):
	""" Returns the degree assortativity coefficient """
	return GraphProperties.degreeAssortativity(G, G.isWeighted())


def properties(G, settings):
	logging.info("[...] calculating properties")

	# size
	n, m = size(G)    # n and m

	logging.info("[...] determining degree distribution")
	# degrees
	degDist = GraphProperties.degreeDistribution(G)
	minDeg, maxDeg, avgDeg = degrees(G)

	plfit = powerLawFit(G, degDist)
	gamma = powerLawExponent(G, degDist)

	# number of isolated nodes
	isolates = degDist[0] if len(degDist) > 0 else None
	satellites = degDist[1] if len(degDist) > 1 else None

	# number of self-loops
	loops = len([(u, v) for (u,v) in G.edges() if (u == v)])

	# density
	dens = density(G)



	# community detection

	ncomPLP, modPLP = None, None
	ncomPLM, modPLM = None, None
	if settings["communities"]:
		logging.info("[...] detecting communities")
		# perform PLM community detection
		logging.info("[...] performing community detection: PLM")
		plm = community.PLM()
		print(plm)
		zetaPLM = plm.run(G)
		ncomPLM = zetaPLM.numberOfSubsets()
		modPLM = community.Modularity().getQuality(zetaPLM, G)

	# degree histogram

	labels, histo = None, None
	if settings["degreeHistogram"]:
		logging.info("[...] preparing degree histogram")
		histo = degDist
		(labels, histo) = compressHistogram(histo, nbins=25)

	# connected components
	nComponents, componentSizes = None, None
	if settings["components"]:
		nComponents, componentSizes = components(G)

	# diameter
	if nComponents == 1:
		dia = Diameter.estimatedDiameterRange(G, error=0.1)
	else:
		dia = None

	# clustering
	avglcc = None
	if settings["clustering"]:
		logging.info("[...] approximating clustering coefficient")
		avglcc = clustering(G)

	# degree assortativity
	logging.info("[...] calculating degree assortativity coefficient")
	assort = degreeAssortativity(G)


	props = {
		 "name": G.getName(),
		 "n": n,
		 "m": m,
		 "minDeg": minDeg,
		 "maxDeg": maxDeg,
		 "avgDeg": avgDeg,
		 "plfit": plfit,
		 "gamma": gamma,
		 "avglcc": avglcc,
		 "nComponents": nComponents,
		 "sizeLargestComponent": max(componentSizes.values()),
		 "dia": dia,
		 "isolates": isolates,
		 "loops": loops,
		 "ncomPLP": ncomPLP,
		 "modPLP": modPLP,
		 "ncomPLM": ncomPLM,
		 "modPLM": modPLM,
		 "dens": dens,
		 "assort": assort,
		 "histo": (labels, histo),
		 }

	return props


def overview(G, settings=collections.defaultdict(lambda: True)):
	"""
	Print an overview of important network properties to the terminal.
	"""
	props = properties(G, settings)
	basicProperties = [
		["nodes (n)", props["n"]],
		["edges (m)", props["m"]],
		["isolated nodes", props["isolates"]],
		["self-loops", props["loops"]],
		["density", "{0:.6f}".format(props["dens"]) if props["dens"] else None],
		["clustering coefficient", "", "{0:.6f}".format(props["avglcc"]) if props["avglcc"] else None],
	]
	degreeProperties = [
		["min. degree", props["minDeg"]],
		["max. degree", props["maxDeg"]],
		["avg. degree", "{0:.6f}".format(props["avgDeg"])],
		["degree power law fit?", "{0}, {1}".format(props["plfit"][0], "{0:.6f}".format(props["plfit"][1]))],
		["degree power law exponent", "{0:.4f}".format(props["gamma"]) if props["plfit"][0] else None],
		["degree assortativity", "{0:.4f}".format(props["assort"])],
	]
	pathStructure = [
		["connected components", props["nComponents"]],
		["size of largest component", props["sizeLargestComponent"]],
		["estimated diameter range", str(props["dia"])],
	]

	communityStructure = [
		["modularity-driven community detection (PLM)", "", ""],
		["", "communities", props["ncomPLM"]],
		["", "modularity", "{0:.6f}".format(props["modPLM"]) if props["modPLM"] else None],
	]

	print()
	print("Network Properties: {0}".format(G.getName()))
	print("==================")
	print("Basic Properties")
	print(tabulate.tabulate(basicProperties))
	print("Node Degree Properties")
	print(tabulate.tabulate(degreeProperties))
	print("Path Structure")
	print(tabulate.tabulate(pathStructure))
	#print("Miscellaneous")
	#print(tabulate.tabulate(miscProperties))
	print("Community Structure")
	print(tabulate.tabulate(communityStructure))
	print("Degree Distribution")
	print("-------------------")
	(labels, histo) = props["histo"]
	if labels and histo:
		termgraph.graph(labels, histo)


def compressHistogram(hist, nbins=20):
	""" Compress a histogram to a number of bins"""
	compressed = [None for i in range(nbins)]
	labels = [None for i in range(nbins)]
	nbinsprev = len(hist)
	binwidth = math.ceil(nbinsprev / nbins)
	for i in range(nbins):
		l = i * binwidth
		u = (i+1) * binwidth
		compressed[i] = sum(hist[l : u])
		labels[i] = "{0}-".format(l)
	return (labels, compressed)



def printDegreeHistogram(G, nbins=25):
	""" Prints a degree histogram as a bar chart to the terminal"""
	hist = GraphProperties.degreeDistribution(G)

	(labels, hist) = compressHistogram(hist, nbins)
	termgraph.graph(labels, hist)
