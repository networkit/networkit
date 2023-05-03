#!/usr/bin/env python3
import unittest
import os
import networkit as nk
import numpy as np
import scipy
import random

class TestCentrality(unittest.TestCase):

	def checkCovers(self, c1, c2):
		if not c1.numberOfElements() == c2.numberOfElements(): return False
		if not c1.numberOfSubsets() == c2. numberOfSubsets(): return False
		for i in range(0,c1.numberOfElements()):
			if not c1.subsetsOf(i) == c2.subsetsOf(i): return False
		return True
	
	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops
		self.LL = nk.readGraph("input/looptest2.gml", nk.Format.GML) #with self-loops sprinkled in
	
	def testApproxBetweenness(self):
		CL = nk.centrality.ApproxBetweenness(self.L, epsilon=0.01, delta=0.1)
		CL.run()
		CLL = nk.centrality.ApproxBetweenness(self.LL, epsilon=0.01, delta=0.1)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))
		for i in range(len(CL.ranking())):
			self.assertAlmostEqual(CL.ranking()[i][1], CLL.ranking()[i][1], delta=0.2*CL.ranking()[i][1])

	def testApproxElectricalCloseness(self):
		for seed in [1, 2, 3]:
			nk.engineering.setSeed(seed, True)
			g = nk.generators.ErdosRenyiGenerator(50, 0.15, False).generate()
			g = nk.components.ConnectedComponents(g).extractLargestConnectedComponent(g, True)
			eps = 0.1
			apx = nk.centrality.ApproxElectricalCloseness(g, eps).run().getDiagonal()

			# Create laplacian matrix
			L = np.zeros((g.numberOfNodes(), g.numberOfNodes()))
			for u in g.iterNodes():
				L[u, u] = g.degree(u)
				for v in g.iterNeighbors(u):
					L[u, v] = -1
					L[v, u] = -1

			pinv = np.linalg.pinv(L).diagonal()
			for u in g.iterNodes():
				self.assertLessEqual(abs(apx[u] - pinv[u]), eps)
	
	def testApproxSpanningEdge(self):
		nk.setSeed(42, False)
		g = nk.generators.ErdosRenyiGenerator(300, 0.1, False).generate()
		g.indexEdges()
		eps = 0.1

		apx = nk.centrality.ApproxSpanningEdge(g, eps)
		apx.run()
		se = nk.centrality.SpanningEdgeCentrality(g, eps)
		se.runParallelApproximation()
		for apxScore, exactScore in zip(apx.scores(), se.scores()):
			self.assertLessEqual(abs(apxScore - exactScore), 2*eps)		

	def testBetweenness(self):
		CL = nk.centrality.Betweenness(self.L)
		CL.run()
		CLL = nk.centrality.Betweenness(self.LL)
		CLL.run()
		self.assertEqual(CL.ranking(), CLL.ranking())

	def testCloseness(self):
		CL = nk.centrality.Closeness(self.L, True, nk.centrality.ClosenessVariant.GENERALIZED)
		CL.run()
		CLL = nk.centrality.Closeness(self.LL, True, nk.centrality.ClosenessVariant.GENERALIZED)
		CLL.run()
		self.assertEqual(CL.ranking(), CLL.ranking())
	
	def testCoreDecomposition(self):
		CL = nk.centrality.CoreDecomposition(self.L)
		CL.run()
		try:
			CLL = nk.centrality.CoreDecomposition(self.LL)
		except RuntimeError:
			import copy
			tmp = copy.deepcopy(self.LL)
			tmp.removeSelfLoops()
			CLL = nk.centrality.CoreDecomposition(tmp)
			CLL.run()
			self.assertTrue(self.checkCovers(CL.getCover(),CLL.getCover()))
	
	def testDegreeCentrality(self):
		g = nk.Graph(8, False, False)

		g.addEdge(0, 2)
		g.addEdge(0, 5)
		g.addEdge(1, 2)
		g.addEdge(2, 3)
		g.addEdge(2, 2)
		g.addEdge(2, 4)
		g.addEdge(3, 5)
		g.addEdge(4, 5)
		g.addEdge(5, 5)
		g.addEdge(5, 6)
		g.addEdge(5, 7)
		g.addEdge(7, 7)

		expected_result = [2.0, 1.0, 4.0, 2.0, 2.0, 5.0, 1.0, 1.0]

		dc = nk.centrality.DegreeCentrality(g).run().scores()

		self.assertListEqual(expected_result, dc)

	def testEigenvectorCentrality(self):
		CL = nk.centrality.EigenvectorCentrality(self.L)
		CL.run()
		CLL = nk.centrality.EigenvectorCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))
	
	def testForest(self):
		nk.engineering.setSeed(42, False)
		eps = 0.05
		g = nk.generators.HyperbolicGenerator(200).generate()
		root = nk.graphtools.augmentGraph(g)

		fc = nk.centrality.ForestCentrality(g, root, eps)
		fc.run()
		apxDiag = fc.getDiagonal()

		A = nk.algebraic.adjacencyMatrix(g, "dense")
		Fmat = scipy.sparse.csgraph.laplacian(A)
		diag = np.linalg.pinv(Fmat).diagonal()

		for apx, exact in zip(apxDiag, diag):
			self.assertLessEqual(abs(apx - exact), eps)

	def testGedWalk(self):
		k, epsilon = 2, 0.05
		gedw = nk.centrality.GedWalk(self.L, k, epsilon)
		gedw.run()
		apxScore, group = gedw.getApproximateScore(), gedw.groupMaxGedWalk()
		self.assertGreaterEqual(apxScore, 0)
		self.assertEqual(len(set(group)), k)
		self.assertAlmostEqual(apxScore, gedw.scoreOfGroup(group), 1)

	def testCentralityGroupClosenessGrowShrink(self):
		g = nk.readGraph('input/MIT8.edgelist', nk.Format.EdgeList, separator='\t', firstNode=0,
				continuous=False, directed=False)
		g = nk.components.ConnectedComponents(g).extractLargestConnectedComponent(g, True)
		k = 5

		nk.engineering.setSeed(42, False)
		for weighted in [False, True]:
			group = set()
			while len(group) < k:
				group.add(nk.graphtools.randomNode(g))

			gc = nk.centrality.GroupClosenessGrowShrink(g, group).run()

			groupMaxCC = gc.groupMaxCloseness()
			self.assertEqual(len(set(groupMaxCC)), k)
			self.assertGreaterEqual(gc.numberOfIterations(), 0)

			for u in groupMaxCC:
				self.assertTrue(g.hasNode(u))

	def testGroupClosenessLocalSearch(self):
		g = nk.readGraph('input/celegans_metabolic.graph', nk.Format.METIS)
		k = 5

		nk.engineering.setSeed(42, False)
		for weighted in [False, True]:
			group = set()
			while len(group) < k:
				group.add(nk.graphtools.randomNode(g))

			gc = nk.centrality.GroupClosenessLocalSearch(g, group).run()

			groupMaxCC = gc.groupMaxCloseness()
			self.assertEqual(len(set(groupMaxCC)), k)

			for u in groupMaxCC:
				self.assertTrue(g.hasNode(u))
	
	def testGroupClosenessLocalSwaps(self):
		k = 5
		g = nk.readGraph('input/MIT8.edgelist', nk.Format.EdgeList, separator='\t', firstNode=0,
				continuous=False, directed=False)
		g = nk.components.ConnectedComponents(g).extractLargestConnectedComponent(g, True)
		for weighted in [False, True]:
			group = set()
			while len(group) < k:
				group.add(nk.graphtools.randomNode(g))
			gc = nk.centrality.GroupClosenessLocalSwaps(g, group).run()

			groupMaxCC = gc.groupMaxCloseness()
			self.assertEqual(len(set(groupMaxCC)), k)
			self.assertGreaterEqual(gc.numberOfSwaps(), 0)

			for u in groupMaxCC:
				self.assertTrue(g.hasNode(u))
	
	def testGroupHarmonicClosenessCentrality(self):
		n, p, k = 50, 0.2, 5
		nk.engineering.setSeed(42, True)
		for directed in [False, True]:
			for weighted in [False, True]:
				g = nk.generators.ErdosRenyiGenerator(n, p, directed).generate()
				if weighted:
					g = nk.graphtools.toWeighted(g)
					g.forEdges(lambda u, v, ew, eid: g.setWeight(u, v, random.random()))

				ghc = nk.centrality.GroupHarmonicCloseness(g, k).run()
				group = ghc.groupMaxHarmonicCloseness()
				self.assertEqual(len(group), k)
				self.assertEqual(len(set(group)), k)
				self.assertGreaterEqual(ghc.scoreOfGroup(g, group), 0)	

	def testKPathCentrality(self):
		CL = nk.centrality.KPathCentrality(self.L)
		CL.run()
		CLL = nk.centrality.KPathCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def testKatzCentrality(self):
		CL = nk.centrality.KatzCentrality(self.L)
		CL.run()
		CLL = nk.centrality.KatzCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))	

	def testPageRank(self):
		CL = nk.centrality.PageRank(self.L)
		CL.run()
		CLL = nk.centrality.PageRank(self.LL)
		CLL.run()

		#test if lists have the same length
		self.assertEqual(len(CL.ranking()), len(CLL.ranking()))
		self.assertEqual(CL.norm, CLL.norm)
		self.assertEqual(CL.maxIterations, CLL.maxIterations)

		norm = nk.centrality.Norm.L2_NORM
		maxIters = 2
		CL.norm = norm
		CL.maxIterations = maxIters
		CLL.norm = norm
		CLL.maxIterations = maxIters

		self.assertEqual(CL.norm, CLL.norm)
		self.assertEqual(CL.maxIterations, CLL.maxIterations)

		CL.run()
		CLL.run()

		self.assertLessEqual(CL.numberOfIterations(), maxIters)
		self.assertLessEqual(CLL.numberOfIterations(), maxIters)
	
	def testRankPerNode(self):
		CL = nk.centrality.PageRank(self.L)
		CL.run()
		CLL = nk.centrality.PageRank(self.LL)
		CLL.run()
		#test if list of pairs and list of ranks have the same length
		self.assertEqual(len(CL.ranking()),len(nk.centrality.rankPerNode(CL.ranking())))
		self.assertEqual(len(CLL.ranking()),len(nk.centrality.rankPerNode(CLL.ranking())))		
	
	def testSquareClusteringCoefficient(self):
		g = nk.Graph()
		g.addNodes(7)
		edges = [(0, 1), (1, 2), (2, 3), (0, 3), (3, 4), (4, 5), (5, 6), (6, 3)]
		[g.addEdge(*edge) for edge in edges]

		expected_result = [1 / 3, 1.0, 1 / 3, 0.2, 1 / 3, 1.0, 1 / 3]
		scores = nk.centrality.LocalSquareClusteringCoefficient(g).run().scores()
		self.assertListEqual(expected_result, scores)

	def testSciPyPageRank(self):
		CL = nk.centrality.SciPyPageRank(self.L)
		CL.run()
		CLL = nk.centrality.SciPyPageRank(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))

	def testSciPyEVZ(self):
		CL = nk.centrality.SciPyEVZ(self.L)
		CL.run()
		CLL = nk.centrality.SciPyEVZ(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))	

	def testRelativeRankErrors(self):
		CL = nk.centrality.Betweenness(self.L)
		CL.run()
		CLL = nk.centrality.Betweenness(self.LL)
		CLL.run()
		self.assertEqual(len(CL.ranking()), len(nk.centrality.relativeRankErrors(CL.ranking(),CLL.ranking())))	
	
	def testTopCloseness(self):
		CC = nk.centrality.Closeness(self.L, True, nk.centrality.ClosenessVariant.GENERALIZED)
		CC.run()
		k = 5
		TC1 = nk.centrality.TopCloseness(self.L, k, True, True)
		TC1.run()
		TC2 = nk.centrality.TopCloseness(self.L, k, True, False)
		TC2.run()
		TC3 = nk.centrality.TopCloseness(self.L, k, False, True)
		TC3.run()
		TC4 = nk.centrality.TopCloseness(self.L, k, False, False)
		TC4.run()
		
		# Test if top nodes and scores lists have the same length
		def testTopKLists(with_trail):
			if not with_trail:
				self.assertEqual(len(TC1.topkNodesList()), k)
				self.assertEqual(len(TC2.topkNodesList()), k)
				self.assertEqual(len(TC3.topkNodesList()), k)
				self.assertEqual(len(TC4.topkNodesList()), k)
			self.assertEqual(len(TC1.topkNodesList(with_trail)), len(TC1.topkScoresList(with_trail)))
			self.assertEqual(len(TC2.topkNodesList(with_trail)), len(TC2.topkScoresList(with_trail)))
			self.assertEqual(len(TC3.topkNodesList(with_trail)), len(TC3.topkScoresList(with_trail)))
			self.assertEqual(len(TC4.topkNodesList(with_trail)), len(TC4.topkScoresList(with_trail)))

		# Test if the ranking is correct
		def testTopKRanking(with_trail):
			def zip_ranking(nodes, scores):
				return [(node, score) for node, score in zip(nodes, scores)]
			length = len(TC1.topkNodesList(with_trail))
			self.assertEqual(CC.ranking()[:length], zip_ranking(TC1.topkNodesList(), TC1.topkScoresList()))
			self.assertEqual(CC.ranking()[:length], zip_ranking(TC2.topkNodesList(), TC2.topkScoresList()))
			self.assertEqual(CC.ranking()[:length], zip_ranking(TC3.topkNodesList(), TC3.topkScoresList()))
			self.assertEqual(CC.ranking()[:length], zip_ranking(TC4.topkNodesList(), TC4.topkScoresList()))

		testTopKLists(False)
		testTopKLists(True)
		testTopKRanking(False)
		testTopKRanking(True)

	def testTopHarmonicCloseness(self):
		CC = nk.centrality.HarmonicCloseness(self.L, False)
		CC.run()
		scores = CC.scores()
		tol = 1e-6
		k = 5

		for useNBbound in [True, False]:
			thc = nk.centrality.TopHarmonicCloseness(self.L, k, useNBbound).run()
			self.assertEqual(len(thc.topkNodesList()), k)
			self.assertEqual(len(thc.topkScoresList()), k)
			self.assertGreaterEqual(len(thc.topkNodesList(True)), k)
			self.assertGreaterEqual(len(thc.topkScoresList(True)), k)

			for node, score in zip(thc.topkNodesList(True), thc.topkScoresList(True)):
				self.assertAlmostEqual(score, scores[node], delta=tol)

			for score_thc, score in zip(thc.topkScoresList(True), sorted(scores, reverse=True)):
				self.assertAlmostEqual(score_thc, score, delta=tol)


if __name__ == "__main__":
	unittest.main()
