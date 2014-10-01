import networkit


# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

class Algo:
	""" runner for an algorithm"""
	def run(self, G):
		raise Exception("Not implemented")

	def loadGraph(self, path):
		G = networkit.readGraph(path, networkit.Format.GML)
		return G

class bConnectedComponents(Algo):
	name = "ConnectedComponents"

	def run(self, G):
		cc = networkit.properties.ConnectedComponents(G)
		cc.run()

class bParallelConnectedComponents(Algo):
	name = "ParallelConnectedComponents"

	def run(self, G):
		cc = networkit.properties.ParallelConnectedComponents(G)
		cc.run()

# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "CoreDecomposition"

	def run(self, G):
		cd = networkit.properties.CoreDecomposition(G)
		cd.run()

# - degree distribution power-law estimation (properties.powerLawExponent)

class bPowerLaw(Algo):
	name = "PowerLaw"

	def run(self, G):
		networkit.properties.powerLawExponent(G)

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(Algo):
	name = "DegreeAssortativity"

	def run(self, G):
		networkit.properties.degreeAssortativity(G)


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
class bBFS(Algo):
	name = "BFS"

	def run(self, G):
		bfs = networkit.graph.BFS(G)
		bfs.run()


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
		d = networkit.properties.Diameter.exactDiameter(G)


class bDiameterEstimate(Algo):
	name = "Diameter"

	def run(self, G):
		d =networkit.properties.Diameter.estimatedDiameterRange(G)

# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

class bClusteringCoefficient(Algo):
	name = "ClusteringCoefficient"

	def run(self, G):
		c = networkit.properties.ClusteringCoefficient.avgLocal(G)

class bApproxClusteringCoefficient(Algo):
	name = "ApproxClusteringCoefficient"

	def run(self, G):
		# TODO: specify number of trials
		c = networkit.properties.ClusteringCoefficient.approxAvgLocal(G)



# - centrality

# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

class bPageRank(Algo):
	name = "PageRank"

	def run(self, G):
		pr = networkit.centrality.PageRank(G)
		pr.run()

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

class bBetweenness(Algo):
	name = "Betweenness"

	def run(self, G):
		bc = networkit.centrality.Betweenness(G)
		bc.run()


class bApproxBetweenness(Algo):
	name = "ApproxBetweenness"

	def run(self, G):
		bc = networkit.centrality.ApproxBetweenness(G, epsilon=0.05, delta=0.1)
		bc.run()
