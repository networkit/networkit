""" This module contains algorithms for the calculation of centrality, i.e. ranking nodes by their structural importance
to the network """


__author__ = "Christian Staudt"

# extension imports
from _NetworKit import Betweenness, PageRank, EigenvectorCentrality, DegreeCentrality, ApproxBetweenness, ApproxBetweenness2, DynBetweenness, DynApproxBetweenness


# local imports
from networkit.algebraic import adjacencyEigenvector, PageRankMatrix, symmetricEigenvectors

# external imports
import math

def ranking(G, algorithm=Betweenness, normalized=False):
	""" Return a ranking of nodes by the specified centrality type"""
	# FIXME: some centrality algorithms take more parameters
	centrality = algorithm(G, normalized)
	centrality.run()
	return centrality.ranking()

def scores(G, algorithm=Betweenness, normalized=False):
	""" Return the centrality scores of nodes using the specified centrality type"""
	centrality = algorithm(G, normalized)
	centrality.run()
	return centrality.scores()


def centralization(G, centralityMeasure):
	"""
	Compute the centralization of a network with respect to some centrality measure.

	The centralization of any network is a measure of how central its most central
	node is in relation to how central all the other nodes are.
	Centralization measures then (a) calculate the sum in differences
	in centrality between the most central node in a network and all other nodes;
	and (b) divide this quantity by the theoretically largest such sum of
	differences in any network of the same size.

	Parameters
	----------
	G : graph
		The graph of which to compute the centrality
	centralityMeasure : instance of Centrality (sub)class
				initialized algorithm object will be run

	Returns
	-------
	double
		centralization score

	"""
	centralityMeasure.run()
	ranking = centralityMeasure.ranking()
	(center, centerScore) = ranking[0]
	maxScore = centralityMeasure.maximum()
	diff1 = sum([(centerScore - c) for (u, c) in ranking])
	diff2 = sum([(maxScore - c) for (u, c) in ranking])
	return diff1 / diff2


class SpectralCentrality:
	"""
	Abstract class to compute the spectral centrality of a graph. This class needs to be supplied with methods
	to generate the correct matrices and do the correct normalization.
	"""
	def __init__(self, G, normalized=False):
		"""
		Constructor.

		Parameters
		----------
		G : graph
		    The graph of which to compute the centrality
		normalized : boolean
					 Whether to normalize the results or not

		"""
		super(SpectralCentrality, self).__init__()

		self.graph = G
		self.normalized = normalized

		self.scoreList = None
		self.rankList = None

	def prepareSpectrum(self):
		""" Method that must be implemented to set the following values:
		self.eigenvectors = list of eigenvectors desired for centrality measure
		self.eigenvalues = list of corresponding eigenvalues
		"""
		raise NotImplemented

	def normFactor(self):
		""" Method that must be implemented to return a correct normalization factor"""
		raise NotImplemented

	def run(self):
		self.prepareSpectrum()

		self.scoreList = None
		self.rankList = None
		self.evz = {}

		if self.normalized:
			normFactor = self.normFactor()
		else:
			normFactor = 1

		for v in self.graph.nodes():
			self.evz[v] = self.eigenvector[v] * normFactor

	def scores(self):
		if self.scoreList is None:
			self.scoreList = [abs(self.evz[v]) for v in self.graph.nodes()] # TODO bad! This depends on iteration order...

		return self.scoreList

	def ranking(self):
		if self.rankList is None:
			self.rankList = [(v, abs(self.evz[v])) for v in self.graph.nodes()]
			self.rankList.sort(key=lambda x: float(x[1]), reverse=True)

		return self.rankList


class SciPyEVZ(SpectralCentrality):
	# TODO: docstring
	def __init__(self, G, normalized=False):
		super(SciPyEVZ, self).__init__(G, normalized=normalized)

	def _length(self, vector):
		square = sum([val * val for val in vector])
		return math.sqrt(square)

	def normFactor(self):
		return 1 / self._length(self.eigenvector)

	def prepareSpectrum(self):
		spectrum = adjacencyEigenvector(self.graph, order=0)
		self.eigenvector = spectrum[1]
		self.eigenvalue = spectrum[0]

class SciPyPageRank(SpectralCentrality):
	# TODO: docstring
	def __init__(self, G, damp=0.95, normalized=False):
		super(SciPyPageRank, self).__init__(G, normalized=normalized)

		self.damp = damp

	def _length(self, vector):
		return sum(vector)

	def normFactor(self):
		return 1 / self._length(self.eigenvector)

	def prepareSpectrum(self):
		prMatrix = PageRankMatrix(self.graph, self.damp)
		spectrum = symmetricEigenvectors(prMatrix, cutoff=0, reverse=False)

		self.eigenvector = spectrum[1][0]
		self.eigenvalue = spectrum[0][0]
