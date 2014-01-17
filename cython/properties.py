# NetworKit native classes and functions
from _NetworKit import GraphProperties, ConnectedComponents

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
	""" Find and analyze detected components """
	logging.info("[...] finding connected components....")
	cc = ConnectedComponents()
	cc.run(G)
	nComponents = cc.numberOfComponents()
	componentSizes = cc.getComponentSizes()
	return (nComponents, componentSizes)

def clustering(G):
	""" Get clustering properties of the graph:
		Returns average local clustering coefficient
	"""
	return GraphProperties.averageLocalClusteringCoefficient(G)


def powerLawExponent(G):
	fit = powerlaw.Fit(degreeDistribution(G))
	return fit.alpha


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
    
	


def properties(nkG, settings):
	nxG = None
	if settings["networkx"]:
		logging.info("[...] converting to NetworX.Graph for some properties....")
		nxG = nxadapter.nk2nx(nkG)

	logging.info("[...] calculating basic properties")
	
	# size
	n, m = size(nkG)    # n and m

	# degrees 
	degDist = GraphProperties.degreeDistribution(nkG) 
	minDeg, maxDeg, avgDeg = degrees(nkG)

	# number of isolated nodes
	isolates = degDist[0] if len(degDist) > 0 else None
	satellites = degDist[1] if len(degDist) > 1 else None

	# number of cliques
	cliques = len(list(nx.find_cliques(nxG))) if nxG else None


	# number of self-loops
	loops = len(nxG.selfloop_edges()) if nxG else None
	
	# density
	dens = density(nkG)



	# community detection

	ncomPLP, modPLP = None, None
	ncomPLM, modPLM = None, None
	if settings["communities"]:
		logging.info("[...] detecting communities")
		# perform PLP and PLM community detection
		logging.debug("performing community detection: PLP")
		# TODO: avoid printout of status bar
		plp = community.PLP()
		zetaPLP = plp.run(nkG)
		ncomPLP = zetaPLP.numberOfClusters()
		modPLP = community.Modularity().getQuality(zetaPLP, nkG)
		logging.info("performing community detection: PLM")
		plm = community.PLM("balanced")
		zetaPLM = plm.run(nkG)
		ncomPLM = zetaPLM.numberOfClusters()
		modPLM = community.Modularity().getQuality(zetaPLM, nkG)

	# degree histogram
	
	labels, histo = None, None
	if settings["degreeDistribution"]:
		logging.info("[...] calculating degree histogram")    
		histo = GraphProperties.degreeDistribution(nkG)
		(labels, histo) = compressHistogram(histo, nbins=25)

	# connected components
	nComponents, componentSizes = None, None
	if settings["components"]:
		logging.info("[...] finding connected components")    
		nComponents, componentSizes = components(nkG)

	# diameter
	dia = None
	if settings["diameter"] and (n < 1000) and (nComponents is 1):
		logging.info("calculating diameter...")
		dia = nx.diameter(nxG)


	# calculate eccentricity
	ecc = None
	if settings["eccentricity"] and (n < 1000) and (nComponents is 1):
		logging.info("calculating eccentricity...")
		eccentricities = nx.eccentricity(nxG)
		ecc = sum(val for val in eccentricities.values()) / n


	# clustering
	avglcc = None
	if settings["clustering"]:
		avglcc = GraphProperties.averageLocalClusteringCoefficient(nkG)

	# degree assortativity
	assort = None
	if settings["assortativity"]:
		assort = nx.degree_assortativity_coefficient(nxG)




	# betweenness centrality
	# TODO: average betweenness centrality?

	props = {
		 "name": nkG.getName(),
		 "n": n,
		 "m": m,
		 "minDeg": minDeg,
		 "maxDeg": maxDeg,
		 "avgDeg": avgDeg,
		 "avglcc": avglcc,
		 "nComponents": nComponents,
		 "sizeLargestComponent": max(componentSizes.values()),
		 "dia": dia,
		 "ecc": ecc,
		 "isolates": isolates,
		 "loops": loops,
		 "ncomPLP": ncomPLP,
		 "modPLP": modPLP,
		 "ncomPLM": ncomPLM,
		 "modPLM": modPLM,
		 "dens": dens,
		 "assort": assort,
		 "cliques": cliques,
		 "histo": (labels, histo),
		 }

	return props


def overview(nkG, settings=collections.defaultdict(lambda: True)):
	"""
	Print an overview of important network properties to the terminal.
	"""
	props = properties(nkG, settings)
	basicProperties = [
		["nodes (n)", props["n"]],
		["edges (m)", props["m"]],
		["min. degree", props["minDeg"]],
		["max. degree", props["maxDeg"]],
		["avg. degree", props["avgDeg"]],
		["isolated nodes", props["isolates"]],
		["self-loops", props["loops"]],
		["density", "{0:.6f}".format(props["dens"]) if props["dens"] else None]
	]
	pathStructure = [
		["connected components", props["nComponents"]],
		["size of largest component", props["sizeLargestComponent"]],
		["diameter", props["dia"]],
		["avg. eccentricity", props["ecc"]],
	]
	
	miscProperties = [
		["degree assortativity", "{0:.6f}".format(props["assort"]) if props["assort"] else None],
		["cliques", props["cliques"]]
	]

	communityStructure = [
		["avg. local clustering coefficient", "", "{0:.6f}".format(props["avglcc"]) if props["avglcc"] else None],
		["PLP community detection", "", ""],
		["", "communities", props["ncomPLP"]],
		["", "modularity", "{0:.6f}".format(props["modPLP"]) if props["modPLP"] else None],
		["PLM community detection", "", ""],
		["", "communities", props["ncomPLM"]],
		["", "modularity", "{0:.6f}".format(props["modPLM"]) if props["modPLM"] else None],
	]

	print()
	print("Network Properties")
	print("==================")
	print("Basic Properties")
	print(tabulate.tabulate(basicProperties))
	print("Path Structure")
	print(tabulate.tabulate(pathStructure))
	print("Miscellaneous")
	print(tabulate.tabulate(miscProperties))
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
	
