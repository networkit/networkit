from _NetworKit import Betweenness, PageRank, EigenvectorCentrality, DegreeCentrality, ApproxBetweenness, ApproxBetweenness2

from algebraic import adjacencyEigenvector

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

class PythonNativeEVZ(object):
	def __init__(self, G, normalized=False):
		super(PythonNativeEVZ, self).__init__()

		self.graph = G
		self.normalized = normalized

		self.scoreList = None
		self.rankList = None

	def prepareSpectrum(self):
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
			self.scoreList = [self.evz[v] for v in self.graph.nodes()] # TODO bad! This depends on iteration order...

		return self.scoreList

	def ranking(self):
		if self.rankList is None:
			self.rankList = [(v, self.evz[v]) for v in self.graph.nodes()]
			self.rankList.sort(key=lambda x: abs(float(x[1])), reverse=True)

		return self.rankList