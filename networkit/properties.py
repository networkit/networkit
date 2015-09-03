# NetworKit native classes and functions
from _NetworKit import ConnectedComponents, ParallelConnectedComponents, StronglyConnectedComponents, ClusteringCoefficient, Diameter, Eccentricity, EffectiveDiameter, Assortativity

# other submodules
from . import community
from . import centrality
from . import termgraph
from . import auxiliary
from . import nxadapter

# other modules
import textwrap
import collections
import math
import logging

try:
	import powerlaw
except ImportError:
	logging.warning("""WARNING: module 'powerlaw' not installed, which is required by some
						functions.""")

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


def degrees(G):
	""" Return min/max/avg degree"""
	avgDeg = GraphProperties.averageDegree(G)
	if G.isDirected():
		minMaxDeg = GraphProperties.minMaxDegreeDirected(G)
	else:
		minMaxDeg = GraphProperties.minMaxDegree(G)
	return (minMaxDeg[0], minMaxDeg[1], avgDeg)

def degreeDistribution(G):
	""" Return the degree distribution of a graph"""
	return GraphProperties.degreeDistribution(G)


def degreeSequence(G):
	""" Return the degree sequence of a graph"""
	return GraphProperties.degreeSequence(G)


def density(G):
	""" Return the density of the graph."""
	(n, m) = G.size()
	loops = G.numberOfSelfLoops()
	m -= loops
	if G.isDirected():
		d = m / (n * (n-1))
	else:
		d = (2 * m) / (n * (n-1))
	return d

# TODO: move to profiling module
def components(G):
	""" Find and analyze detected components.
		Returns the number of components and the sizes
		of each component. For more details use the
		ConnectedComponents class.
	"""
	logging.info("[...] finding connected components....")
	if G.isDirected():
		cc = StronglyConnectedComponents(G)
	else:
		cc = ConnectedComponents(G)
	cc.run()
	components = cc.getPartition()
	nComponents = components.numberOfSubsets()
	componentSizes = components.subsetSizeMap()
	return (nComponents, componentSizes)

def numberOfComponents(G):
	""" Find and number of components """
	logging.info("[...] finding connected components....")
	cc = ConnectedComponents(G)
	cc.run()
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
	if G.numberOfNodes() < 100:
		return ClusteringCoefficient().avgLocal(G)
	else:
		nSamples = math.ceil(math.log(10) / (error**2)) # fixed confidence of 90%
		logging.info("taking {0} samples".format(nSamples))
		return ClusteringCoefficient().approxAvgLocal(G, nSamples)


def degreePowerLaw(G, dd=None):
	""" Check if a power law is a good fit for the degree distribution.

	Returns
	-------
	answer: bool
		whether a power law is a good fit
	R : double
		goodness of the fit, i.e.
		the loglikelihood ratio between the two candidate distributions. This number will be positive if the data is more likely in the first distribution, and negative if the data is more likely in the 		  	 	second distribution. The exponential distribution is the absolute minimum alternative candidate for evaluating the heavy- tailedness of the distribution. The reason is definitional: the typical quantitative definition of a ”heavy- tail” is that it is not exponentially bounded. Thus if a power law is not a better fit than an exponential distribution (as in the above example) there is scarce ground for considering the distribution to be heavy-tailed at all, let alone a power law.
	gamma : double
		the degree power law exponent
	"""
	if not dd:
		dd = degreeSequence(G)
	fit = powerlaw.Fit(dd)
	R, p = fit.distribution_compare("power_law", "exponential", normalized_ratio=True)
	gamma = fit.alpha
	return ((R > 0), R, gamma)



def degreeAssortativity(G):
	""" Returns the degree assortativity coefficient """
	if G.isDirected():
		results = []
		results.append(GraphProperties.degreeAssortativityDirected(G,0,0))
		results.append(GraphProperties.degreeAssortativityDirected(G,1,1))
		return results
	else:
		return GraphProperties.degreeAssortativity(G, G.isWeighted())


def degeneracy(G):
	""" degeneracy of an undirected graph is defined as the largest k for which
	the graph has a non-empty k-core. Degeneracy is only implemented for graphs without
	self-loops."""
	if G.numberOfSelfLoops() > 0:
		raise NotImplementedError("Call Graph.removeSelfLoops() first.")
	coreDec = centrality.CoreDecomposition(G)
	coreDec.run()
	return coreDec.maxCoreNumber()
