"""
This module provides graph generators that produce synthetic networks according to various models.
"""

__author__ = "Christian Staudt"

# extension imports
from _NetworKit import Graph, BarabasiAlbertGenerator, PubWebGenerator, ErdosRenyiGenerator, ClusteredRandomGraphGenerator, DorogovtsevMendesGenerator, DynamicPubWebGenerator, DynamicPathGenerator, ChungLuGenerator, HyperbolicGenerator, DynamicHyperbolicGenerator, HavelHakimiGenerator, DynamicDorogovtsevMendesGenerator, RmatGenerator, DynamicForestFireGenerator, RegularRingLatticeGenerator, WattsStrogatzGenerator, PowerlawDegreeSequence, EdgeSwitchingMarkovChainGenerator, EdgeSwitchingMarkovChainGenerator as ConfigurationModelGenerator, LFRGenerator


class StarGraphGenerator:
	"""
	Creates a graph with star topology, i.e. one
	central node connected to all other nodes.
	"""

	def __init__(self, n):
		self.n = n

	def generate(self):
		G = Graph(self.n)
		for u in G.nodes():
			if u is 0:
				pass
			else:
				G.addEdge(u, 0)
		return G


class MultiscaleGenerator:
	"""
	TODO:
	"""

	def __init__(self, O, withADWeights=True):
		self.O = O				# original graph
		self.Gc = []			# hierarchy of coarse graphs
		self.Gf = []			# hierarchy of fine graphs
		self.matching = []
		self.up = []			# mapping: fine node -> coarse node
		self.down = []			# mapping: coarse node -> fine node
		self.nodeWeights = []

	def _buildCoarseSequence(self, maxLevel):
		for i in range(maxLevel):
			if i is 0:
				self.Gc.append(self.O)
				self.up.append({})
				self.down.append({})
				self.nodeWeights.append([1 for v in range(self.Gc[i].upperNodeIdBound())])
			else:
				if self.withADWeights:
					ad = distance.AlgebraicDistance(self.Gc[i], withEdgeScores=True).preprocess()
					matcher = matching.PathGrowingMatcher(self.Gc[i], weightsFromADScores(ad.getEdgeScores(), self.Gc[i]))
				else:
					matcher = matching.PathGrowingMatcher(self.Gc[i])
				self.matching.append(matcher.run().getMatching())
				coarsening = coarsening.MatchingCoarsening(self.Gc[i], self.matching[i], noSelfLoops=True)
				coarsening.run()
				self.Gc.append(coarsening.getCoarseGraph())
				self.Gc[i].indexEdges()
				self.up[i], self.down[i] = coarsening.getFineToCoarseNodeMapping(), coarsening.getCoarseToFineNodeMapping()
				if (self.Gc[i].size() == self.Gc[h-1].size()):
					print("stopping at level ", h)
					break
				else:
					pass
					#nodeWeights = [0 for v in range(G.upperNodeIdBound())]
					#self.Gc[i].forNodes(lambda v: nodeWeights[v] += sum())

	def generate(self):
		return None

	@classmethod
	def fit(cls, G):
		return cls(G)
