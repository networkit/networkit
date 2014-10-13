try:
	import igraph
except ImportError as ex:
	print("igraph not available")

import random


class Algo:

	frameworkPrefix = "ig:"

	""" runner for an algorithm"""
	def run(self, G):
		raise Exception("Not implemented")

	def loadGraph(self, path):
		G = igraph.read(path, format="gml")  # Check if the format is actually edgelist
		return G


class bConnectedComponents(Algo):
	name = "ig:ConnectedComponents"

	def run(self, G):
		comp = G.components()
		return len(comp.sizes())

# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "ig:CoreDecomposition"

	def run(self, G):
		igraph.Graph.coreness(G)


# - degree distribution power-law estimation (properties.powerLawExponent)

class bPowerLaw(Algo):
	name = "ig:PowerLaw"

	def run(self, G):
		fit = igraph.statistics.power_law_fit(G.degree())
		return fit.alpha

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(Algo):
	name = "ig:DegreeAssortativity"

	def run(self, G):
		return igraph.Graph.assortativity_degree(G)


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)

class bBFS(Algo):
	name = "ig:BFS"

	def run(self, G):
		s = random.randint(0, G.vcount())
		G.bfs(s)


# - community detection (community.PLM, community.PLP)
class bCommunityDetectionLM(Algo):
	name = "ig:CommunityDetectionLM"

	def run(self, G):
		igraph.Graph.community_multilevel(G)


class bCommunityDetectionLP(Algo):
	name = "ig:CommunityDetectionLP"

	def run(self, G):
		igraph.Graph.community_label_propagation(G)


# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)


class bDiameter(Algo):
	name = "ig:Diameter"

	def run(self, G):
		return igraph.Graph.diameter(G)


class bDiameterEstimate(Algo):
	name = "ig:Diameter"

	def run(self, G):
		raise NotImplementedError("TODO")  # Can't find relevant method


# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)
class bClusteringCoefficient(Algo):
	name = "ig:ClusteringCoefficient"

	def run(self, G):
		return igraph.Graph.transitivity_local_undirected(G)


class bApproxClusteringCoefficient(Algo):
	name = "ig:ApproxClusteringCoefficient"

	def run(self, G):
		# TODO: specify number of trials
		raise NotImplementedError("TODO")  # Can't find relevant method

# - centrality


# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)
class bPageRank(Algo):
	name = "ig:PageRank"

	def run(self, G):
		pr = igraph.Graph.pagerank(G)	# this leads to a segfault -seems buggy

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)
class bBetweenness(Algo):
	name = "ig:Betweenness"

	def run(self, G):
		b = igraph.Graph.betweenness(G)


class bApproxBetweenness(Algo):
	name = "ig:ApproxBetweenness"

	def run(self, G):
		raise NotImplementedError("TODO")  # Can't find relevant method
