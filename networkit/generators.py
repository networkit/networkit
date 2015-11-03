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
