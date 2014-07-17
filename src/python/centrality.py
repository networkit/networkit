from _NetworKit import Betweenness, PageRank, EigenvectorCentrality, DegreeCentrality, ApproxBetweenness, ApproxBetweenness2

from algebraic import adjacencyEigenvector, PageRankMatrix, symmetricEigenvectors

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

class SpectralCentrality(object):		
	def __init__(self, G, normalized=False):
		super(SpectralCentrality, self).__init__()

		self.graph = G
		self.normalized = normalized

		self.scoreList = None
		self.rankList = None

	def prepareSpectrum(self):
		raise NotImplemented
		spectrum = adjacencyEigenvector(self.graph, order=0)
		self.eigenvector = spectrum[1]
		self.eigenvalue = spectrum[0]

	def length(self, vector):
		square = sum([val * val for val in vector])
		return math.sqrt(square)

	def run(self):
		self.prepareSpectrum()

		self.scoreList = None
		self.rankList = None
		self.evz = {}

		if self.normalized:
			normFactor = 1 / self.length(self.eigenvector)
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


class PythonNativeEVZ(SpectralCentrality):
	def __init__(self, G, normalized=False):
		super(PythonNativeEVZ, self).__init__(G, normalized=normalized)

	def prepareSpectrum(self):
		spectrum = adjacencyEigenvector(self.graph, order=0)
		self.eigenvector = spectrum[1]
		self.eigenvalue = spectrum[0]

class PythonNativePageRank(SpectralCentrality):
	def __init__(self, G, damp=0.95, normalized=False):
		super(PythonNativePageRank, self).__init__(G, normalized=normalized)

		self.damp = damp

	def prepareSpectrum(self):
		prMatrix = PageRankMatrix(self.graph, self.damp)
		spectrum = symmetricEigenvectors(prMatrix, cutoff=0, reverse=False)

		self.eigenvector = spectrum[1][0]
		self.eigenvalue = spectrum[0][0]