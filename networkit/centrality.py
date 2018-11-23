""" This module contains algorithms for the calculation of centrality, i.e. ranking nodes by their structural importance
to the network """


__author__ = "Christian Staudt"
__credits__ = ["Christian Staudt", "Elisabetta Bergamini", "Henning Meyerhenke", "Marc Nemes", "Maximilian Vogel"]

# extension imports
# TODO: (+) ApproxCloseness

from _NetworKit import Betweenness, PageRank, EigenvectorCentrality, DegreeCentrality, ApproxBetweenness,\
ApproxBetweenness2, EstimateBetweenness, DynApproxBetweenness, Closeness, HarmonicCloseness, KPathCentrality, CoreDecomposition,\
KatzCentrality, LocalClusteringCoefficient, ApproxCloseness, LocalPartitionCoverage, Sfigality, SpanningEdgeCentrality,\
PermanenceCentrality, TopCloseness, TopHarmonicCloseness, DynTopHarmonicCloseness, DynBetweenness,\
GroupDegree, GroupCloseness, DynBetweennessOneNode, LaplacianCentrality, ApproxGroupBetweenness, DynKatzCentrality,\
KadabraBetweenness
from _NetworKit import _ClosenessVariant as ClosenessVariant

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



def rankPerNode(ranking):
	"""
	Parameters
	----------
 	ranking: ordered list of tuples (node, score)

	Returns
	_______
	for each node (sorted by node ID), the ranking of the node

	"""
	n_nodes = len(ranking)
	ranking_id = [0]*n_nodes
	for index, pair in enumerate(ranking):
		ranking_id[pair[0]] = index
	#we assign to all nodes the ranking of the first node with the same score
	for index, pair in enumerate(ranking):
			if index == 0:
				continue
			if pair[1] == ranking[index-1][1]:
				prev_node = ranking[index-1][0]
				ranking_id[pair[0]] = ranking_id[prev_node]
	return ranking_id


def relativeRankErrors(rx, ry):
	"""
	Let $r_x(u)$ be the rank of node $u$ in ranking $x$.
	The relative rank error of node $u$ is defined as
		$$r_x(u) / r_y(u)$$


	Parameters
	----------
	rx : list
		ranking - ordered list of tuples (node, score)

	ry:  list
		ranking - ordered list of tuples (node, score)

	Returns
	_______
	list of rank errors ordered by node ID

	"""
	diff = []
	n = len(rx)
	if not(n == len(ry)):
		return diff
	rnode_x = rankPerNode(rx)
	rnode_y = rankPerNode(ry)
	for i in range(n):
		diff.append((rnode_x[i]+1)/(rnode_y[i]+1))
	return diff


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
		self.evz = {}

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
		return self

	def scores(self):
		if self.scoreList is None:
			self.scoreList = [v for k,v in self.evz.items()]

		return self.scoreList

	def ranking(self):
		if self.rankList is None:
			self.rankList = sorted(self.evz.items(),key=lambda x: float(x[1]), reverse=True)
		return self.rankList


class SciPyEVZ(SpectralCentrality):
	"""
	Compute Eigenvector centrality using algebraic meh

	Parameters
	----------
	G : graph
		The graph of which to compute the centrality
	normalized : boolean
				 Whether to normalize the results or not

	"""
	def __init__(self, G, normalized=False):
		if G.isDirected():
			raise NotImplementedError("Not implemented for directed graphs; use centrality.EigenvectorCentrality instead")
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
