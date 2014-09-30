try:
	import networkx
except ImportError:
	error("networkx not available")


class Algo:
	""" runner for an algorithm"""
	def run(self, G):
		raise Exception("Not implemented")

	def loadGraph(self, path):
		G = networkx.read_edgelist(path)
		return G

	def numberOfEdges(self, G):
		return G.numberOfEdges()


# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

class bConnectedComponents(Algo):
	name = "ConnectedComponents"

	def run(self, G):
		cc = networkx.connected_components(G)

# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "CoreDecomposition"

	def run(self, G):
		cn = networkx.core_number(G)

# - degree distribution power-law estimation (properties.powerLawExponent)

# not available

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(Algo):
	name = "DegreeAssortativity"

	def run(self, G):
		ac = networkx.degree_assortativity_coefficient(G)


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
class bBFS(Algo):
	name = "BFS"

	def run(self, G):
		networkx.bfs_predecessors(G)


# - community detection (community.PLM, community.PLP)
# not available

# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)

class bDiameter(Algo):
	name = "Diameter"

	def run(self, G):
		d = networkx.diameter(G)


# approximate diameter not available

# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

class bClusteringCoefficient(Algo):
	name = "ClusteringCoefficient"

	def run(self, G):
		c = properties.ClusteringCoefficient.avgLocal(G)

class bApproxClusteringCoefficient(Algo):
	name = "ApproxClusteringCoefficient"

	def run(self, G):
		# TODO: number of trials
		c = networkx.average_clustering(G, trials)



# - centrality

# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

class bPageRank(Algo):
	name = "PageRank"

	def run(self, G):
		# TODO: check parameters
		pr = networkx.pagerank(G)

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

class bBetweenness(Algo):
	name = "Betweenness"

	def run(self, G):
		networkx.betweenness_centrality(G)


# approximation not available

#-------------------------
