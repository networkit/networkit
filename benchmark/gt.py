import graph_tool
import random


# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

framework = "(gt)"


class Algo:

	framework = framework

	""" runner for an algorithm"""
	def run(self, G):
		raise Exception("Not implemented")

	def loadGraph(self, path):
		G = graph_tool.load_graph(path, fmt="gml")
		return G

class bConnectedComponents(Algo):
	name = "ConnectedComponents" + framework

	def run(self, G):
		cc = graph_tool.topology.label_components(G)



# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "CoreDecomposition" + framework

	def run(self, G):
		cd = graph_tool.topology.kcore_decomposition(G)

# - degree distribution power-law estimation (properties.powerLawExponent)

# class bPowerLaw(Algo) – not found in GT

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(Algo):
	name = "DegreeAssortativity" + framework

	def run(self, G):
		graph_tool.correlations.scalar_assortativity(G, deg="total")


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
class bBFS(Algo):
	name = "BFS" + framework

	def run(self, G):
		v = random.randint(0, G.num_vertices()-1)
		graph_tool.topology.shortest_distance(G, v)

# - community detection (community.PLM, community.PLP)

#class bCommunityDetectionLM(Algo): – not implemented in graph-tool

#class bCommunityDetectionLP(Algo): – not implemented in graph-tool

# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)

class bDiameter(Algo):
	name = "Diameter" + framework

	def run(self, G):
		d = graph_tool.stats.distance_histogram(G)


class bDiameterEstimate(Algo):
	name = "DiameterEstimate" + framework

	def run(self, G):
		d = graph_tool.topology.pseudo_diameter(G)

# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

class bClusteringCoefficient(Algo):
	name = "ClusteringCoefficient" + framework

	def run(self, G):
		c = graph_tool.clustering.local_clustering(G)

# class bApproxClusteringCoefficient(Algo): – not implemented in graph-tool



# - centrality

# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

class bPageRank(Algo):
	name = "PageRank" + framework

	def run(self, G):
		pr = graph_tool.centrality.pagerank(G)

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)

class bEigenvector(Algo):
	name = "Eigenvector" + framework

	def run(self, G):
		pr = graph_tool.eigenvector(G)

# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

class bBetweenness(Algo):
	name = "Betweenness" + framework

	def run(self, G):
		graph_tool.centrality.betweenness(G)


#class bApproxBetweenness(Algo):
#	name = "ApproxBetweenness"

#	def run(self, G):
#		graph_tool.betweenness(G)
