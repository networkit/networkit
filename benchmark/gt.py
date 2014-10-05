import graph_tool.all as graphtool


# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

class Algo:
	""" runner for an algorithm"""
	def run(self, G):
		raise Exception("Not implemented")

	def loadGraph(self, path):
		G = graphtool.load_graph(path, fmt="gml")
		return G

class bConnectedComponents(Algo):
	name = "ConnectedComponents"

	def run(self, G):
		cc = graphtool.label_components(G)



# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "CoreDecomposition"

	def run(self, G):
		cd = graphtool.kcore_decomposition(G)

# - degree distribution power-law estimation (properties.powerLawExponent)

# class bPowerLaw(Algo) – not found in GT

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(Algo):
	name = "DegreeAssortativity"

	def run(self, G):
		graphtool.assortativity(G, deg="total")


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
class bBFS(Algo):
	name = "BFS"

	def run(self, G):
		gt.shortest_distance(G) # Johnson’s algorithm seems to be used https://graph-tool.skewed.de/static/doc/topology.html#graph_tool.topology.shortest_distance


# - community detection (community.PLM, community.PLP)

class bCommunityDetectionLM(Algo):
	name = "CommunityDetectionLM"

	def run(self, G):
		plm = networkit.community.PLM(G)
		plm.run()

class bCommunityDetectionLP(Algo):
	name = "CommunityDetectionLP"

	def run(self, G):
		plm = networkit.community.PLP(G)
		plm.run()

# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)

class bDiameter(Algo):
	name = "Diameter"

	def run(self, G):
		d = graphtool.distance_histogram(G)


class bDiameterEstimate(Algo):
	name = "Diameter"

	def run(self, G):
		d = graphtool.pseudo_diameter(G)

# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

class bClusteringCoefficient(Algo):
	name = "ClusteringCoefficient"

	def run(self, G):
		c = graphtool.local_clustering(G)

# class bApproxClusteringCoefficient(Algo): – not implemented in graph-tool



# - centrality

# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

class bPageRank(Algo):
	name = "PageRank"

	def run(self, G):
		pr = graphtool.pagerank(G)

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)

class bEigenvector(Algo):
	name = "Eigenvector"

	def run(self, G):
		pr = graphtool.eigenvector(G)

# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

class bBetweenness(Algo):
	name = "Betweenness"

	def run(self, G):
		graphtool.betweenness(G)


#class bApproxBetweenness(Algo):
#	name = "ApproxBetweenness"

#	def run(self, G):
#		graphtool.betweenness(G)
