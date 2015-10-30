import networkit

from util import *
import base

framework = "(nk)"

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


# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "CoreDecomposition" + framework

	def run(self, G):
		cd = networkit.centrality.CoreDecomposition(G)
		cd.run()

class bCoreDecompositionSeq(Algo):
	name = "CoreDecompositionSeq" + framework

	def run(self, G):
		cd = networkit.centrality.CoreDecomposition(G, enforceBucketQueueAlgorithm=True)
		cd.run()


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
		plm = networkit.community.PLM(G, turbo=True)
		plm.run()

class bCommunityDetectionLP(Algo):
	name = "CommunityDetectionLP" + framework

	def run(self, G):
		plm = networkit.community.PLP(G)
		plm.run()

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
		networkit.centrality.LocalClusteringCoefficient(G).run()


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

class bEigenvectorCentrality(Algo):
	name = "EigenvectorCentrality" + framework

	def run(self, G):
		evc = networkit.centrality.EigenvectorCentrality(G, tol=1e-06)
		evc.run()

class bKatzCentrality(Algo):
	name = "KatzCentrality" + framework

	def run(self, G):
		kc = networkit.centrality.KatzCentrality(G, tol=1e-06)
		kc.run()


class bDegreeAssortativity(Algo):
	name = "DegreeAssortativity" + framework

	def run(self, G):
		networkit.correlation.Assortativity(G, networkit.centrality.DegreeCentrality(G).run().scores()).run()


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

class bBetweenness(Algo):
	name = "Betweenness" + framework

	def run(self, G):
		bc = networkit.centrality.Betweenness(G)
		bc.run()

class bBetweennessSeq(Algo):
	name = "BetweennessSeq" + framework

	def run(self, G):
		mt = networkit.getMaxNumberOfThreads()
		networkit.setNumberOfThreads(1)
		bc = networkit.centrality.Betweenness(G)
		bc.run()
		networkit.setNumberOfThreads(mt)


class bApproxBetweenness(Algo):
	name = "BetweennessApprox" + framework

	def run(self, G):
		bc = networkit.centrality.ApproxBetweenness2(G, nSamples=42)
		bc.run()


class bApproxBetweennessSeq(Algo):
	name = "BetweennessApproxSeq" + framework

	def run(self, G):
		mt = networkit.getMaxNumberOfThreads()
		networkit.setNumberOfThreads(1)
		bc = networkit.centrality.ApproxBetweenness2(G, nSamples=42)
		bc.run()
		networkit.setNumberOfThreads(mt)
