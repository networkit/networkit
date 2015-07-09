try:
	import igraph
except ImportError as ex:
	print("igraph not available")

import random

framework = " (ig)"


class Algo:

	framework = framework

	""" runner for an algorithm"""
	def run(self, G):
		raise Exception("Not implemented")

	def loadGraph(self, path):
		G = igraph.read(path, format="gml")  # Check if the format is actually edgelist
		return G


class bConnectedComponents(Algo):
	name = "ConnectedComponents" + framework

	def run(self, G):
		comp = G.components()
		return len(comp.sizes())

# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "CoreDecomposition" + framework

	def run(self, G):
		igraph.Graph.coreness(G)


# - degree distribution power-law estimation (properties.powerLawExponent)

class bPowerLaw(Algo):
	name = "PowerLaw" + framework

	def run(self, G):
		fit = igraph.statistics.power_law_fit(G.degree())
		return fit.alpha

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(Algo):
	name = "DegreeAssortativity" + framework

	def run(self, G):
		return igraph.Graph.assortativity_degree(G)


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)

class bBFS(Algo):
	name = "BFS" + framework

	def run(self, G):
		s = random.randint(0, G.vcount())
		G.bfs(s)


# - community detection (community.PLM, community.PLP)
class bCommunityDetectionLM(Algo):
	name = "CommunityDetectionLM" + framework

	def run(self, G):
		igraph.Graph.community_multilevel(G)


class bCommunityDetectionLP(Algo):
	name = "CommunityDetectionLP" + framework

	def run(self, G):
		igraph.Graph.community_label_propagation(G)


# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)


class bDiameter(Algo):
	name = "Diameter" + framework

	def run(self, G):
		return igraph.Graph.diameter(G)



# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)
class bClusteringCoefficient(Algo):
	name = "ClusteringCoefficient" + framework

	def run(self, G):
		igraph.Graph.transitivity_local_undirected(G)



# - centrality


# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)
class bPageRank(Algo):
	name = "PageRank" + framework

	def run(self, G):
		pr = igraph.Graph.pagerank(G)	# this leads to a segfault -seems buggy

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)
class bBetweenness(Algo):
	name = "Betweenness" + framework

	def run(self, G):
		b = igraph.Graph.betweenness(G)
