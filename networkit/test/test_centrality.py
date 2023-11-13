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

	def testApproxBetweennessNumberOfSamples(self):
		CL = nk.centrality.ApproxBetweenness(self.L, epsilon=0.01, delta=0.1)
		CL.run()
		self.assertEqual(CL.numberOfSamples(), 63026)

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
		self.assertEqual(len(CLL.edgeScores()), 0)
		self.assertAlmostEqual(CLL.centralization(), 1.22, 2)
		self.assertEqual(CLL.maximum(), 28.0)

	def testCloseness(self):
		CL = nk.centrality.Closeness(self.L, True, nk.centrality.ClosenessVariant.GENERALIZED)
		CL.run()
		CLL = nk.centrality.Closeness(self.LL, True, nk.centrality.ClosenessVariant.GENERALIZED)
		CLL.run()
		self.assertEqual(CL.ranking(), CLL.ranking())

	def testClosenessApprox(self):
		#expecting same results from exact algorithm and approx with 50 samples
		apr = nk.centrality.ApproxCloseness(self.L, 50, True, nk.centrality.ClosenessVariant.GENERALIZED)
		apr.run()
		ex=nk.centrality.Closeness(self.L, True, nk.centrality.ClosenessVariant.GENERALIZED)
		ex.run()
		self.assertListEqual(apr.ranking(), ex.ranking())	
		self.assertListEqual(apr.getSquareErrorEstimates(), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])	
	
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
   
	def testComplexPathsAllNodes(self):
		CP = nk.centrality.ComplexPaths(self.L, 3)
		CP.run()

		res = CP.getPLci() 
		self.assertEqual(len(res), self.L.numberOfNodes())     
  
	def testComplexPathsSingleNode(self):
		G = nk.Graph(8)
		G.addEdge(0, 2)
		G.addEdge(1, 2)
		G.addEdge(2, 3)
		G.addEdge(2, 4)
		G.addEdge(3, 5)
		G.addEdge(4, 5)
		G.addEdge(4, 7)
		G.addEdge(5, 6)
		G.addEdge(6, 7)
  
		CP = nk.centrality.ComplexPaths(G, 2, "singleNode", 0)
		CP.run()

		res = CP.getAdopters()
		self.assertListEqual(res, [0,2])
       
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

	def testDynBetweennessRun(self):
		#results of dynBetweenness against Betweenness should be equal
		dyn = nk.centrality.DynBetweenness(self.L)
		ex = nk.centrality.Betweenness(self.L)
		dyn.run()
		ex.run()
		self.assertListEqual(dyn.ranking(), ex.ranking())
		self.assertListEqual(dyn.scores(), ex.scores())
		self.assertEqual(dyn.score(0), ex.score(0))

	def testDynBetweennessUpdate(self):
		G = nk.Graph(4)
		G.addEdge(0,1)
		dyn = nk.centrality.DynBetweenness(G)
		dyn.run()
		up1 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 0, 3, 4.0)	
		dyn.update(up1)
		self.assertListEqual(dyn.ranking(), [(0, 2.0), (1, 0.0), (2, 0.0), (3, 0.0)])

	def testDynBetweennessUpdateBatch(self):
		G = nk.Graph(4)
		G.addEdge(0,1)
		dyn = nk.centrality.DynBetweenness(G)
		dyn.run()	
		up2 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 2, 3, 1.0)	
		up3 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 1, 3, 1.0)	
		dyn.updateBatch([up2,up3])
		self.assertListEqual(dyn.ranking(), [(1, 2.0), (0, 0.0), (2, 0.0), (3, 0.0)])

	def testDynBetweennessApproxRun(self):
		dyn = nk.centrality.DynApproxBetweenness(self.L)
		dyn.run()
		#due to float inaccuracies we only compare the length of results instead of values
		self.assertEqual(len(dyn.ranking()), 9)
		self.assertEqual(len(dyn.scores()), 9)
		self.assertEqual(dyn.score(0), 0.0)
		self.assertEqual(dyn.getNumberOfSamples(), 63026)

	"""def testDynBetweennessApproxUpdate(self):
		G = nk.Graph(4)
		G.addEdge(0,1)
		dyn = nk.centrality.DynApproxBetweenness(G).run()
		up1 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 0, 3, 4.0)	
		dyn.update(up1)
		self.assertListEqual(dyn.ranking(), [(0, 0.0), (1, 0.0), (2, 0.0), (3, 0.0)])		

	def testDynBetweennessApproxUpdateBatch(self):
		G = nk.Graph(4)
		G.addEdge(0,1)
		dyn = nk.centrality.DynApproxBetweenness(G).run()		
		up2 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 2, 3, 1.0)	
		up3 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 1, 3, 1.0)	
		dyn.updateBatch([up2,up3])
		self.assertListEqual(dyn.ranking(), [(0, 0.0), (1, 0.0), (2, 0.0), (3, 0.0)])"""		

	#calling run currently results in segfault
	"""def testDynBetweennessOneNodeUpdate(self):
		G = nk.Graph(5)
		G.addEdge(0,1)
		dynOne = nk.centrality.DynBetweennessOneNode(G, 4)
		#dynOne.run()
		up1 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 0, 2, 1.0)	
		dynOne.update(up1)
		self.assertEqual(dynOne.getDistance(0,3), 1.0)
		self.assertEqual(dynOne.getSigma(0,3), 1.0)
		self.assertEqual(dynOne.getSigmax(0,3), 1.0)
		self.assertEqual(dynOne.getbcx(), 4.0)
		
	def testDynBetweennessOneNodeUpdates(self):
		G = nk.Graph(4)
		G.addEdge(0,1)
		G.addEdge(2,3)
		dynOne = nk.centrality.DynBetweennessOneNode(G, 0)
		dynOne.run()		
		up2 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 0, 2, 1.0)	
		up3 = nk.dynamics.GraphEvent(nk.dynamics.GraphEventType.EDGE_ADDITION, 0, 3, 1.0)	
		dynOne.updateBatch([up2,up3])
		self.assertEqual(dynOne.getDistance(0,3), 0.0)
		self.assertEqual(dynOne.getSigma(0,3), 1.0)
		self.assertEqual(dynOne.getSigmax(0,3), 1.0)
		self.assertEqual(dynOne.getbcx(), 4.0)	"""
	
	def testDynKatzCentrality(self):
		CL = nk.centrality.DynKatzCentrality(self.L, 9)
		CL.run()
		#due to float inaccuracies we only compare the length of results instead of values
		self.assertEqual(len(CL.ranking()), 9)	
		self.assertEqual(CL.top(), 4)
		self.assertAlmostEqual(CL.bound(0), 0.62, 2)
		self.assertTrue(CL.areDistinguished(0,8))

	def testDynTopHarmonicCloseness(self):
		k = 5
		
		thc = nk.centrality.DynTopHarmonicCloseness(self.L, k, False).run()
		hc = nk.centrality.HarmonicCloseness(self.L, False).run()
		score = thc.ranking()
		refScore = hc.ranking()
		
		#asserting equal first k values, respecing float inaccuracies
		for i in range(k):
			self.assertAlmostEqual(score[i][1], refScore[i][1], delta = 0.5)
		self.assertListEqual(thc.topkNodesList(), [4, 1, 6, 2, 3])

		topK = [5.833, 5.416, 5.333, 4.916, 4.916]
		for i in range(k):
			self.assertAlmostEqual(thc.topkScoresList()[i], topK[i], delta=0.5)

	def testEigenvectorCentrality(self):
		CL = nk.centrality.EigenvectorCentrality(self.L)
		CL.run()
		CLL = nk.centrality.EigenvectorCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))

	def testEstimateBetweeness(self):
		CL = nk.centrality.EstimateBetweenness(self.L, 50)
		CL.run()
		self.assertEqual(len(CL.ranking()), 9)		
	
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
	
	def testGroupBetweennessApprox(self):
		gc = nk.centrality.ApproxGroupBetweenness(self.L, 3, 0.5)
		gc.run()
		self.assertEqual(len(gc.groupMaxBetweenness()), 3)

	def testGroupBetweennessApproxScoreOfGroup(self):
		gc = nk.centrality.ApproxGroupBetweenness(self.L, 3, 0.5)	
		self.assertEqual(gc.scoreOfGroup([0,1,2]), 14.0)
	
	def testGroupCloseness(self):
		gc = nk.centrality.GroupCloseness(self.L)
		gc.run()
		self.assertListEqual(gc.groupMaxCloseness(), [4])

	def testGroupClosenessFarness(self):
		gc = nk.centrality.GroupCloseness(self.L)
		self.assertEqual(gc.computeFarness([0,1,2]), 14.0)

	def testGroupClosenessScoreOfGroup(self):	
		gc = nk.centrality.GroupCloseness(self.L)
		self.assertAlmostEqual(gc.scoreOfGroup([0,1,2]), 0.43, 2)
	
	def testGroupDegree(self):
		gd = nk.centrality.GroupDegree(self.L)
		gd.run()
		self.assertEqual(gd.getScore(), 5)
		self.assertListEqual(gd.groupMaxDegree(), [4])

	def testGroupDegreeScoreOfGroup(self):	
		gd = nk.centrality.GroupDegree(self.L)
		self.assertEqual(gd.scoreOfGroup([0,1,2]), 5)
	
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
			self.assertGreaterEqual(gc.numberOfIterations(), 0)
	
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

	def testKadabraBetweenness(self):
		KB = nk.centrality.KadabraBetweenness(self.L, err=0.01, delta=0.1, k=3)
		KB.run()
		#check ranking() returns correct rankings
		for i in range(len(KB.ranking())):
			self.assertEqual(KB.ranking()[i][0], KB.topkNodesList()[i])
			self.assertAlmostEqual(KB.ranking()[i][1], KB.topkScoresList()[i], delta=0.1)
		
		#check scores() returns correct scores
		scores = [0.0, 0.397, 0.0, 0.0, 0.896, 0.0, 0.844, 0.395, 0.0]
		for i in range(self.L.numberOfNodes()):
			self.assertAlmostEqual(KB.scores()[i], scores[i], delta=0.1)

		#check helper functions, which are non deterministic
		self.assertGreaterEqual(KB.getNumberOfIterations(), 10000)
		self.assertGreaterEqual(KB.getOmega(), 10000)
		
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

	def testLaplacianCentrality(self):
		LC = nk.centrality.LaplacianCentrality(self.L)
		LC.run()
		LCC = nk.centrality.LaplacianCentrality(self.LL)
		LCC.run()
		self.assertEqual(len(LC.ranking()),len(LCC.ranking()))

	
	def testLocalClustetingCoefficient(self):
		LCC = nk.centrality.LocalClusteringCoefficient(self.L)
		LCC.run()
		LCCC = nk.centrality.LaplacianCentrality(self.LL)
		LCCC.run()
		self.assertEqual(len(LCC.ranking()),len(LCCC.ranking()))		

	def testLocalPartitionCoverage(self):
		part=nk.structures.Partition(9)
		part.addToSubset(0, 1)
		part.addToSubset(0, 2)
		part.addToSubset(0, 3)
		part.addToSubset(1, 4)
		part.addToSubset(1, 5)
		CL = nk.centrality.LocalPartitionCoverage(self.L, part)
		CL.run()
	
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
	
	def testPermanenceCentrality(self):
		part=nk.structures.Partition(9)
		part.addToSubset(0, 1)
		part.addToSubset(0, 2)
		part.addToSubset(0, 3)
		part.addToSubset(1, 4)
		part.addToSubset(1, 5)
		PC = nk.centrality.PermanenceCentrality(self.L, part)
		PC.run()
		self.assertEqual(PC.getIntraClustering(0), 0.0)
		self.assertEqual(PC.getPermanence(0), -1.0)
	
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
