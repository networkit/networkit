import networkit

from util import *
import base

framework = " (nk)"

# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

class Algo(base.Algo):
	""" runner for an algorithm"""

	framework = framework

	def loadGraph(self, path):
		with Timer() as t:
			G = networkit.readGraph(path, networkit.Format.GML)
		debug("reading {path} took {t.elapsed} s".format(**locals()))
		return G

class bConnectedComponents(Algo):
	name = "ConnectedComponents" + framework

	def run(self, G):
		cc = networkit.properties.ConnectedComponents(G)
		cc.run()
		return cc.numberOfComponents()

class bParallelConnectedComponents(Algo):
	name = "ParallelConnectedComponents" + framework

	def run(self, G):
		cc = networkit.properties.ParallelConnectedComponents(G)
		cc.run()
		return cc.numberOfComponents()

# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "CoreDecomposition" + framework

	def run(self, G):
		cd = networkit.properties.CoreDecomposition(G)
		cd.run()

# - degree distribution power-law estimation (properties.powerLawExponent)

class bPowerLaw(Algo):
	name = "PowerLaw" + framework

	def run(self, G):
		return networkit.properties.degreePowerLaw(G)

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(Algo):
	name = "DegreeAssortativity" + framework

	def run(self, G):
		return networkit.properties.degreeAssortativity(G)


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
class bBFS(Algo):
	name = "BFS" + framework

	def run(self, G):
		bfs = networkit.graph.BFS(G, G.randomNode(), storePaths=False)
		bfs.run()


# - community detection (community.PLM, community.PLP)

class bCommunityDetectionLM(Algo):
	name = "CommunityDetectionLM" + framework

	def run(self, G):
		plm = networkit.community.PLM(turbo=True)
		plm.run(G)

class bCommunityDetectionLP(Algo):
	name = "CommunityDetectionLP" + framework

	def run(self, G):
		plm = networkit.community.PLP()
		plm.run(G)

# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)

class bDiameter(Algo):
	name = "Diameter" + framework

	def run(self, G):
		return networkit.properties.Diameter.exactDiameter(G)


class bDiameterEstimate(Algo):
	name = "DiameterEstimate" + framework

	def run(self, G):
		return networkit.properties.Diameter.estimatedDiameterRange(G)

# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

class bClusteringCoefficient(Algo):
	name = "ClusteringCoefficient" + framework

	def run(self, G):
		c = networkit.properties.ClusteringCoefficient.avgLocal(G)
		return c

class bApproxClusteringCoefficient(Algo):
	name = "ClusteringCoefficientApprox" + framework

	def run(self, G):
		# TODO: specify number of trials
		c = networkit.properties.ClusteringCoefficient.approxAvgLocal(G, trials=1000)
		return c



# - centrality

# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

class bPageRank(Algo):
	name = "PageRank" + framework

	def run(self, G):
		pr = networkit.centrality.PageRank(G, damp=0.85, tol=1e-06)
		pr.run()

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

class bBetweenness(Algo):
	name = "Betweenness" + framework

	def run(self, G):
		bc = networkit.centrality.Betweenness(G)
		bc.run()


class bApproxBetweenness(Algo):
	name = "BetweennessApprox" + framework

	def run(self, G):
		bc = networkit.centrality.ApproxBetweenness(G, epsilon=0.1, delta=0.1)
		bc.run()
