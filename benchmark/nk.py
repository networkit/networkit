import networkit

from util import *
import base

framework = "networkit"

# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

class Algo(base.Algo):
	""" runner for an algorithm"""

	framework = framework

	def loadGraph(self, path, graphFormat=networkit.Format.GML):
		with Timer() as t:
			G = networkit.readGraph(path, graphFormat)
		debug("reading {path} took {t.elapsed} s".format(**locals()))
		return G

class bConnectedComponents(Algo):
	name = "ConnectedComponents"

	def run(self, G):
		cc = networkit.components.ConnectedComponents(G)
		cc.run()
		return cc.numberOfComponents()


# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(Algo):
	name = "CoreDecomposition"

	def run(self, G):
		cd = networkit.centrality.CoreDecomposition(G)
		cd.run()

class bCoreDecompositionSeq(Algo):
	name = "CoreDecompositionSeq"

	def run(self, G):
		cd = networkit.centrality.CoreDecomposition(G, enforceBucketQueueAlgorithm=True)
		cd.run()


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
class bBFS(Algo):
	name = "BFS"

	def run(self, G):
		bfs = networkit.graph.BFS(G, G.randomNode(), storePaths=False)
		bfs.run()


# - community detection (community.PLM, community.PLP)

class bCommunityDetectionLM(Algo):
	name = "CommunityDetectionLM"

	def run(self, G):
		plm = networkit.community.PLM(G, turbo=True)
		plm.run()

class bCommunityDetectionLP(Algo):
	name = "CommunityDetectionLP"

	def run(self, G):
		plm = networkit.community.PLP(G)
		plm.run()


# distance module

class bDiameter(Algo):
	name = "Diameter"

	def run(self, G):
		return networkit.distance.Diameter(G, networkit.distance.DiameterAlgo.exact).run()


class bDiameterEstimate(Algo):
	name = "DiameterEstimate"

	def run(self, G):
		return networkit.distance.Diameter(G, networkit.distance.DiameterAlgo.estimatedRange, error=0.1).run()


class bEffectiveDiameter(Algo):
	name = "EffectiveDiameter"

	def run(self, G):
		return networkit.distance.EffectiveDiameter(G).run().getEffectiveDiameter()


class bEffectiveDiameterApproximation(Algo):
	name = "EffectiveDiameterApproximation"

	def run(self, G):
		return networkit.distance.EffectiveDiameterApproximation(G).run().getEffectiveDiameter()


class bApproxHopPlot(Algo):
	name = "ApproxHopPlot"

	def run(self, G):
		return networkit.distance.ApproxHopPlot(G).run().getHopPlot()


class bNeighborhoodFunction(Algo):
	name = "NeighborhoodFunction"

	def run(self, G):
		return networkit.distance.NeighborhoodFunction(G).run().getNeighborhoodFunction()


class bApproxNeighborhoodFunction(Algo):
	name = "ApproxNeighborhoodFunction"

	def run(self, G):
		return networkit.distance.ApproxNeighborhoodFunction(G).run().getNeighborhoodFunction()



# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

class bClusteringCoefficient(Algo):
	name = "ClusteringCoefficient"

	def run(self, G):
		networkit.centrality.LocalClusteringCoefficient(G).run()


class bApproxClusteringCoefficient(Algo):
	name = "ClusteringCoefficientApprox"

	def run(self, G):
		# TODO: specify number of trials
		c = networkit.properties.ClusteringCoefficient.approxAvgLocal(G, trials=1000)
		return c



# - centrality

# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

class bPageRank(Algo):
	name = "PageRank"

	def run(self, G):
		pr = networkit.centrality.PageRank(G, damp=0.85, tol=1e-06)
		pr.run()

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)

class bEigenvectorCentrality(Algo):
	name = "EigenvectorCentrality"

	def run(self, G):
		evc = networkit.centrality.EigenvectorCentrality(G, tol=1e-06)
		evc.run()

class bKatzCentrality(Algo):
	name = "KatzCentrality"

	def run(self, G):
		kc = networkit.centrality.KatzCentrality(G, tol=1e-06)
		kc.run()


class bDegreeAssortativity(Algo):
	name = "DegreeAssortativity"

	def run(self, G):
		networkit.correlation.Assortativity(G, networkit.centrality.DegreeCentrality(G).run().scores()).run()


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness)

class bBetweenness(Algo):
	name = "Betweenness"

	def run(self, G):
		bc = networkit.centrality.Betweenness(G)
		bc.run()

class bBetweennessSeq(Algo):
	name = "BetweennessSeq"

	def run(self, G):
		mt = networkit.getMaxNumberOfThreads()
		networkit.setNumberOfThreads(1)
		bc = networkit.centrality.Betweenness(G)
		bc.run()
		networkit.setNumberOfThreads(mt)


class bApproxBetweenness(Algo):
	name = "BetweennessApprox"

	def run(self, G):
		bc = networkit.centrality.EstimateBetweenness(G, nSamples=42)
		bc.run()


class bApproxBetweennessSeq(Algo):
	name = "BetweennessApproxSeq"

	def run(self, G):
		mt = networkit.getMaxNumberOfThreads()
		networkit.setNumberOfThreads(1)
		bc = networkit.centrality.EstimateBetweenness(G, nSamples=42)
		bc.run()
		networkit.setNumberOfThreads(mt)
