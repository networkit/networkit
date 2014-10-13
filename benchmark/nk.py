import networkit

from util import *
import base

# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

class Algo(base.Algo):
	""" runner for an algorithm"""

	frameworkPrefix = "nk:"

	def loadGraph(self, path):
		with Timer() as t:
			G = networkit.readGraph(path, networkit.Format.GML)
		debug("reading {path} took {t.elapsed} s".format(**locals()))
		return G

class bConnectedComponents(Algo):
	name = "nk:ConnectedComponents"

	def run(self, G):
		cc = networkit.properties.ConnectedComponents(G)
		cc.run()
		return cc.numberOfComponents()

class bParallelConnectedComponents(Algo):
	name = "nk:ParallelConnectedComponents"

	def run(self, G):
		cc = networkit.properties.ParallelConnectedComponents(G)
		cc.run()
		return cc.numberOfComponents()

# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "nk:CoreDecomposition"

	def run(self, G):
		cd = networkit.properties.CoreDecomposition(G)
		cd.run()

# - degree distribution power-law estimation (properties.powerLawExponent)

class bPowerLaw(Algo):
	name = "nk:PowerLaw"

	def run(self, G):
		return networkit.properties.powerLawExponent(G)

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(Algo):
	name = "nk:DegreeAssortativity"

	def run(self, G):
		return networkit.properties.degreeAssortativity(G)


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
class bBFS(Algo):
	name = "nk:BFS"

	def run(self, G):
		bfs = networkit.graph.BFS(G, G.randomNode())
		bfs.run()


# - community detection (community.PLM, community.PLP)

class bCommunityDetectionLM(Algo):
	name = "nk:CommunityDetectionLM"

	def run(self, G):
		plm = networkit.community.PLM()
		plm.run(G)

class bCommunityDetectionLP(Algo):
	name = "nk:CommunityDetectionLP"

	def run(self, G):
		plm = networkit.community.PLP()
		plm.run(G)

# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)

class bDiameter(Algo):
	name = "nk:Diameter"

	def run(self, G):
		return networkit.properties.Diameter.exactDiameter(G)


class bDiameterEstimate(Algo):
	name = "nk:DiameterEstimate"

	def run(self, G):
		return networkit.properties.Diameter.estimatedDiameterRange(G)

# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

class bClusteringCoefficient(Algo):
	name = "nk:ClusteringCoefficient"

	def run(self, G):
		c = networkit.properties.ClusteringCoefficient.avgLocal(G)
		return c

class bApproxClusteringCoefficient(Algo):
	name = "nk:ApproxClusteringCoefficient"

	def run(self, G):
		# TODO: specify number of trials
		c = networkit.properties.ClusteringCoefficient.approxAvgLocal(G, trials=1000)
		return c



# - centrality

# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

class bPageRank(Algo):
	name = "nk:PageRank"

	def run(self, G):
		pr = networkit.centrality.PageRank(G, damp=0.85, tol=1e-06)
		pr.run()

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

class bBetweenness(Algo):
	name = "nk:Betweenness"

	def run(self, G):
		bc = networkit.centrality.Betweenness(G)
		bc.run()


class bApproxBetweenness(Algo):
	name = "nk:ApproxBetweenness"

	def run(self, G):
		bc = networkit.centrality.ApproxBetweenness(G, epsilon=0.2, delta=0.1)
		bc.run()
