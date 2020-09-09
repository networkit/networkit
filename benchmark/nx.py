try:
	import networkx
except ImportError:
	error("networkx not available")

import random

framework = " (nx)"

class Algo:
	framework = framework

	""" runner for an algorithm"""
	def run(self, G):
		raise Exception("Not implemented")

	def loadGraph(self, path):
		G = networkx.read_gml(path)
		return G

	def preprocessGraph(self, G):
		return G	# do nothing

	def numberOfEdges(self, G):
		return G.numberOfEdges()


# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

class bConnectedComponents(Algo):
	name = "ConnectedComponents" + framework

	def run(self, G):
		cc = networkx.number_connected_components(G)
		return cc

# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "CoreDecomposition" + framework

	def run(self, G):
		cn = networkx.core_number(G)
		#return cn[0]

# - degree distribution power-law estimation (properties.powerLawExponent)

# not available

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(Algo):
	name = "DegreeAssortativity" + framework

	def run(self, G):
		ac = networkx.degree_assortativity_coefficient(G)
		return ac


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
class bBFS(Algo):
	name = "BFS" + framework

	def run(self, G):
		s = random.randint(0, G.number_of_nodes())
		networkx.bfs_predecessors(G, s)


# - community detection (community.PLM, community.PLP)
# not available

# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)

class bDiameter(Algo):
	name = "Diameter" + framework

	def run(self, G):
		d = networkx.diameter(G)
		return d


# approximate diameter not available

# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

class bClusteringCoefficient(Algo):
	name = "ClusteringCoefficient" + framework

	def run(self, G):
		c = networkx.average_clustering(G)
		return c

class bApproxClusteringCoefficient(Algo):
	name = "ApproxClusteringCoefficient" + framework

	def run(self, G):
		c = networkx.average_clustering(G, trials=1000)
		return C



# - centrality

# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

class bPageRank(Algo):
	name = "PageRank" + framework

	def run(self, G):
		pr = networkx.pagerank(G, alpha=0.85, tol=1e-06)

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness)

class bBetweenness(Algo):
	name = "Betweenness" + framework

	def run(self, G):
		networkx.betweenness_centrality(G)


# approximation not available

#-------------------------
