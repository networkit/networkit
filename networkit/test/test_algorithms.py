#!/usr/bin/env python3
import numpy as np
import os
import random
import unittest

import networkit as nk

class TestSelfLoops(unittest.TestCase):

	def checkCovers(self, c1, c2):
		if not c1.numberOfElements() == c2.numberOfElements(): return False
		if not c1.numberOfSubsets() == c2. numberOfSubsets(): return False
		for i in range(0,c1.numberOfElements()):
			if not c1.subsetsOf(i) == c2.subsetsOf(i): return False
		return True

	def setUp(self):
		# toggle the comment/uncomment to test on small or large test cases
		#self.L = nk.readGraph("PGPgiantcompo.graph", nk.Format.METIS) #without self-loops
		#self.LL = nk.readGraph("PGPConnectedCompoLoops.gml", nk.Format.GML) #with self-loops sprinkled in
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops
		self.LL = nk.readGraph("input/looptest2.gml", nk.Format.GML) #with self-loops sprinkled in

	def testCentralityBetweenness(self):
		CL = nk.centrality.Betweenness(self.L)
		CL.run()
		CLL = nk.centrality.Betweenness(self.LL)
		CLL.run()
		self.assertEqual(CL.ranking(), CLL.ranking())

	def testCentralityApproxBetweenness(self):
		CL = nk.centrality.ApproxBetweenness(self.L, epsilon=0.01, delta=0.1)
		CL.run()
		CLL = nk.centrality.ApproxBetweenness(self.LL, epsilon=0.01, delta=0.1)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))
		for i in range(len(CL.ranking())):
			self.assertAlmostEqual(CL.ranking()[i][1], CLL.ranking()[i][1], delta=0.2*CL.ranking()[i][1])


	def testCentralityCloseness(self):
		CL = nk.centrality.Closeness(self.L, True, nk.centrality.ClosenessVariant.Generalized)
		CL.run()
		CLL = nk.centrality.Closeness(self.LL, True, nk.centrality.ClosenessVariant.Generalized)
		CLL.run()
		self.assertEqual(CL.ranking(), CLL.ranking())


	def testCentralityTopCloseness(self):
		CC = nk.centrality.Closeness(self.L, True, nk.centrality.ClosenessVariant.Generalized)
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

	def testCentralityTopHarmonicCloseness(self):
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

	def testCentralityCoreDecomposition(self):
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


	def testCentralityEigenvectorCentrality(self):
		CL = nk.centrality.EigenvectorCentrality(self.L)
		CL.run()
		CLL = nk.centrality.EigenvectorCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def testCentralityKPathCentrality(self):
		CL = nk.centrality.KPathCentrality(self.L)
		CL.run()
		CLL = nk.centrality.KPathCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def testCentralityKatzCentrality(self):
		CL = nk.centrality.KatzCentrality(self.L)
		CL.run()
		CLL = nk.centrality.KatzCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def testCentralityPageRank(self):
		CL = nk.centrality.PageRank(self.L)
		CL.run()
		CLL = nk.centrality.PageRank(self.LL)
		CLL.run()

		#test if lists have the same length
		self.assertEqual(len(CL.ranking()), len(CLL.ranking()))
		self.assertEqual(CL.norm, CLL.norm)
		self.assertEqual(CL.maxIterations, CLL.maxIterations)

		norm = nk.centrality.Norm.l2norm
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

	def testCentralityRankPerNode(self):
		CL = nk.centrality.PageRank(self.L)
		CL.run()
		CLL = nk.centrality.PageRank(self.LL)
		CLL.run()
		#test if list of pairs and list of ranks have the same length
		self.assertEqual(len(CL.ranking()),len(nk.centrality.rankPerNode(CL.ranking())))
		self.assertEqual(len(CLL.ranking()),len(nk.centrality.rankPerNode(CLL.ranking())))

	def testGedWalkCentrality(self):
		k, epsilon = 2, 0.05
		gedw = nk.centrality.GedWalk(self.L, k, epsilon)
		gedw.run()
		apxScore, group = gedw.getApproximateScore(), gedw.groupMaxGedWalk()
		self.assertGreaterEqual(apxScore, 0)
		self.assertEqual(len(set(group)), k)
		self.assertAlmostEqual(apxScore, gedw.scoreOfGroup(group), 1)

	def test_centrality_SciPyPageRank(self):
		CL = nk.centrality.SciPyPageRank(self.L)
		CL.run()
		CLL = nk.centrality.SciPyPageRank(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def testCentralitySciPyEVZ(self):
		CL = nk.centrality.SciPyEVZ(self.L)
		CL.run()
		CLL = nk.centrality.SciPyEVZ(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))

	def testCentralityRelativeRankErrors(self):
		CL = nk.centrality.Betweenness(self.L)
		CL.run()
		CLL = nk.centrality.Betweenness(self.LL)
		CLL.run()
		self.assertEqual(len(CL.ranking()), len(nk.centrality.relativeRankErrors(CL.ranking(),CLL.ranking())))

	def testCentralityApproxSpanningEdge(self):
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

	def testCommunityPLM(self):
		PLML = nk.community.PLM(self.L)
		PLMLL = nk.community.PLM(self.LL)
		PLML.run()
		PLMLL.run()
		PLMLP = PLML.getPartition()
		PLMLLP = PLMLL.getPartition()
		self.assertIsNot(PLMLP.getSubsetIds(), None)
		self.assertIsNot(PLMLLP.getSubsetIds(), None)
		# test if partitions add up to original set
		reconstructedSet = []
		for i in PLMLP.getSubsetIds():
			for j in PLMLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.L.iterNodes()), set(reconstructedSet))
		reconstructedSet = []
		for i in PLMLLP.getSubsetIds():
			for j in PLMLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.LL.iterNodes()), set(reconstructedSet))

	def testCommunityPLP(self):
		PLPL = nk.community.PLP(self.L)
		PLPLL = nk.community.PLP(self.LL)
		PLPL.run()
		PLPLL.run()
		PLPLP = PLPL.getPartition()
		PLPLLP = PLPLL.getPartition()
		self.assertIsNot(PLPLP.getSubsetIds(), None)
		self.assertIsNot(PLPLLP.getSubsetIds(), None)
		# test if partitions add up to original set
		reconstructedSet = []
		for i in PLPLP.getSubsetIds():
			for j in PLPLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.L.iterNodes()), set(reconstructedSet))
		reconstructedSet = []
		for i in PLPLLP.getSubsetIds():
			for j in PLPLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.LL.iterNodes()), set(reconstructedSet))


	def testCommunityCutClustering(self):
		CL = nk.community.CutClustering(self.L, 0.2)
		CLL = nk.community.CutClustering(self.LL, 0.2)
		CL.run()
		CLL.run()
		CLP = CL.getPartition()
		CLLP = CLL.getPartition()
		self.assertIsNot(CLP.getSubsetIds(), None)
		self.assertIsNot(CLLP.getSubsetIds(), None)
		# test if partitions add up to original set
		reconstructedSet = []
		for i in CLP.getSubsetIds():
			for j in CLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.L.iterNodes()), set(reconstructedSet))
		reconstructedSet = []
		for i in CLLP.getSubsetIds():
			for j in CLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.LL.iterNodes()), set(reconstructedSet))


	def testCommunityGraphClusteringTools(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		GCT = nk.community.GraphClusteringTools()
		self.assertIsInstance(GCT.equalClustering(PLMLLP,PLPLLP, self.LL), bool)
		self.assertIsInstance(GCT.getImbalance(PLMLLP), float)
		self.assertIsInstance(GCT.isOneClustering(self.LL, PLPLLP), bool)
		self.assertIsInstance(GCT.isProperClustering(self.LL, PLMLLP), bool)
		self.assertIsInstance(GCT.isSingletonClustering(self.LL, PLPLLP), bool)


	def testCommunityGraphStructuralRandMeasure(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		GSRM = nk.community.GraphStructuralRandMeasure()
		self.assertAlmostEqual(GSRM.getDissimilarity(self.LL, PLMLLP, PLPLLP),0.5, delta=0.5 )


	def testCommunityHubdominance(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		HD = nk.community.HubDominance()
		self.assertIsInstance(HD.getQuality(PLMLLP, self.LL),float )


	def testCommunityJaccardMeasure(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		JM = nk.community.JaccardMeasure()
		self.assertIsInstance(JM.getDissimilarity(self.LL, PLMLLP, PLPLLP),float)


	def testCommunityLPDegreeOrdered(self):
		CL = nk.community.LPDegreeOrdered(self.L)
		CLL = nk.community.LPDegreeOrdered(self.LL)
		CL.run()
		CLL.run()
		CLP = CL.getPartition()
		CLLP = CLL.getPartition()
		self.assertIsNot(CLP.getSubsetIds(), None)
		self.assertIsNot(CLLP.getSubsetIds(), None)
		# test if partitions add up to original set
		reconstructedSet = []
		for i in CLP.getSubsetIds():
			for j in CLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.L.iterNodes()), set(reconstructedSet))
		reconstructedSet = []
		for i in CLLP.getSubsetIds():
			for j in CLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.LL.iterNodes()), set(reconstructedSet))

	def testCommunityModularity(self):
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		Mod = nk.community.Modularity()
		self.assertAlmostEqual(Mod.getQuality(PLPLLP, self.LL),0.25, delta=0.75)


	def testCommunityNMIDistance(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		NMI = nk.community.NMIDistance()
		self.assertIsInstance(NMI.getDissimilarity(self.LL, PLMLLP, PLPLLP),float)


	def testCommunityNodeStructuralRandMeasure(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		NSRM = nk.community.NodeStructuralRandMeasure()
		self.assertAlmostEqual(NSRM.getDissimilarity(self.LL, PLMLLP, PLPLLP),0.5, delta=0.5 )


	def testCommunityCommunityGraph(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		CG = nk.community.communityGraph(self.LL, PLMLLP)


	def testCommunityEvaluateCommunityDetection(self):
		PLMLL = nk.community.PLM(self.LL)
		nk.community.evalCommunityDetection(PLMLL, self.LL)


	def testCommunityKCoreCommunityDetection(self):
		with self.assertRaises(RuntimeError) as cm:
			kCCD = nk.community.kCoreCommunityDetection(self.LL, 1, inspect=False)


	def testFlowEdmondsKarp(self):
		self.L.indexEdges()
		self.LL.indexEdges()
		r1 = nk.graphtools.randomNode(self.L)
		r2 = nk.graphtools.randomNode(self.L)
		while r1 is r2:
			r2 = nk.graphtools.randomNode(self.L)
		EKL = nk.flow.EdmondsKarp(self.L, r1, r2)
		EKLL = nk.flow.EdmondsKarp(self.LL, r1, r2)
		EKL.run()
		EKLL.run()


	def testGlobalsClusteringCoefficient(self):
		CL = nk.globals.ClusteringCoefficient()
		CL.exactGlobal(self.L)
		CL.exactGlobal(self.LL)
		CL.approxGlobal(self.L, 5)
		CL.approxGlobal(self.LL, 5)
		CL.approxAvgLocal(self.L, 5)
		CL.approxAvgLocal(self.LL, 5)
		CL.avgLocal(self.L)
		with self.assertRaises(RuntimeError) as cm:
			CL.avgLocal(self.LL)
		CL.sequentialAvgLocal(self.L)
		CL.sequentialAvgLocal(self.LL)


	def testComponentsConnectedComponents(self):
		CC = nk.components.ConnectedComponents(self.LL)
		CC.run()
		CC.componentOfNode(1)
		CC.getComponentSizes()
		CC.getPartition()
		CC.numberOfComponents()

	def testExtractLargestConnectedComponent(self):
		G = nk.Graph(10)
		for i in range(3):
			G.addEdge(i, i+1)

		for i in range(4, 9):
			G.addEdge(i, i+1)

		G1 = nk.components.ConnectedComponents.extractLargestConnectedComponent(G, True)
		self.assertEqual(G1.numberOfNodes(), 6)
		self.assertEqual(G1.numberOfEdges(), 5)

		G2 = nk.components.ConnectedComponents.extractLargestConnectedComponent(G, False)
		for i in range(G.numberOfNodes()):
			self.assertEqual(G2.hasNode(i), (4 <= i <= 9))

	def test_components_StronglyConnectedComponents(self):
		g = nk.readGraph("input/MIT8.edgelist",
				nk.Format.EdgeList, separator='\t', firstNode=0, continuous=False, directed=True)
		scc = nk.components.StronglyConnectedComponents(g)
		scc.run()
		self.assertNotEqual(scc.componentOfNode(0), None)
		nComponents = scc.numberOfComponents()
		compSizes = scc.getComponentSizes()
		self.assertEqual(nComponents, len(compSizes))

		comps = scc.getComponents()
		for idx, size in compSizes.items():
			self.assertEqual(len(comps[idx]), size)

		_=scc.getPartition()

	def testComponentsBiconnectedComponents(self):
		bcc = nk.components.BiconnectedComponents(self.LL)
		bcc.run()

		for component in bcc.getComponents():
			G1 = nk.graphtools.subgraphFromNodes(self.LL, component)
			def testNode(v):
				G2 = nk.Graph(G1)
				G2.removeNode(v)
				cc = nk.components.ConnectedComponents(G2)
				cc.run()
				self.assertEqual(cc.numberOfComponents(), 1)
			G1.forNodes(testNode)


	def testDistanceDiameter(self):
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.EstimatedRange, error = 0.1)
		D.run()
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.EstimatedSamples, nSamples = 5)
		D.run()
		D = nk.distance.Diameter(self.LL, nk.distance.DiameterAlgo.Exact)
		D.run()


	def testDistanceEccentricity(self):
		E = nk.distance.Eccentricity()
		E.getValue(self.LL, 0)


	def testDistanceEffectiveDiameter(self):
		algo = nk.distance.EffectiveDiameter(self.L)
		algo.run()
		algo = nk.distance.EffectiveDiameter(self.LL)
		algo.run()


	def testDistanceApproxEffectiveDiameter(self):
		algo = nk.distance.EffectiveDiameterApproximation(self.L)
		algo.run()
		algo = nk.distance.EffectiveDiameterApproximation(self.LL)
		algo.run()


	def testDistanceApproxHopPlot(self):
		algo = nk.distance.HopPlotApproximation(self.L)
		algo.run()
		algo = nk.distance.HopPlotApproximation(self.LL)
		algo.run()


	def testDistanceNeighborhoodFunction(self):
		algo = nk.distance.NeighborhoodFunction(self.L)
		algo.run()
		algo = nk.distance.NeighborhoodFunction(self.LL)
		algo.run()


	def testDistanceApproxNeighborhoodFunction(self):
		algo = nk.distance.NeighborhoodFunctionApproximation(self.L)
		algo.run()
		algo = nk.distance.NeighborhoodFunctionApproximation(self.LL)
		algo.run()

	def testDistanceAStar(self):
		# Builds a mesh graph with the given number of rows and columns
		def buildMesh(rows, cols):
			G = nk.Graph(rows * cols, False, False)
			for i in range(rows):
				for j in range(cols):
					if j < cols - 1:
						G.addEdge(i * cols + j, i * cols + j + 1)
					if i < rows - 1:
						G.addEdge(i * cols + j, (i + 1) * cols + j)
			return G

		# Test the AStar algorithm on a mesh with the given number of rows and columns
		def testMesh(rows, cols):
			G = buildMesh(rows, cols)

			# Test A* on the given source-target pair
			def testPair(s, t):

				# Some distance heuristics:

				# Always returns 0, A* degenerates to Dijkstra
				def zeroDist(u):
					return 0

				# Returns the exact distance from u to the target
				def exactDist(u):
					rowU = int(u / cols)
					colU = int(u % cols)
					rowT = int(t / cols)
					colT = int(t % cols)
					return abs(rowU - rowT) + abs(colU - colT)

				# Returns the eucledian distance from u to the target
				def eucledianDist(u):
					rowT = int(t / cols)
					colT = int(t % cols)
					rowDiff = abs(int(u / cols) - rowT)
					colDiff = abs(int(u % cols) - rowT)
					return (rowDiff**2 + colDiff**2)**.5

				# Use BFS as ground truth
				bfs = nk.distance.BFS(G, s, True, False, t).run()

				# Test A* on all the heuristics
				for heu in [zeroDist, exactDist, eucledianDist]:
					heuristics = [heu(u) for u in range(G.numberOfNodes())]
					astar = nk.distance.AStar(G, heuristics, s, t, True)
					astar.run()

					# Test distance of target
					self.assertEqual(astar.getDistance(), bfs.distance(t))

					# Test path
					path = astar.getPath()
					self.assertEqual(len(path), len(bfs.getPath(t)) - 2)
					if len(path) == 0:
						continue
					for i in range(len(path) - 1):
						self.assertTrue(G.hasEdge(path[i], path[i + 1]))

			# Iterate over all possible source-target pairs
			G.forNodePairs(testPair)

		# Test some meshes
		testMesh(10, 10)
		testMesh(21, 5)
		testMesh(9, 18)
		testMesh(7, 1)

	def testApproxElectricalClosenes(self):
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

	def testDistanceSPSP(self):
		nk.engineering.setSeed(1, True)
		random.seed(1)
		for directed in [True, False]:
			for weighted in [True, False]:
				g = nk.generators.ErdosRenyiGenerator(100, 0.15, directed).generate()
				if weighted:
					g = nk.graphtools.toWeighted(g)
					g.forEdges(lambda u, v, ew, eid: g.setWeight(u, v, random.random()))
				for nSources in [1, 10, 50]:
					sources = set()
					while len(sources) < nSources:
						sources.add(nk.graphtools.randomNode(g))

					spsp = nk.distance.SPSP(g, list(sources))
					spsp.run()
					spsp.getDistances()

if __name__ == "__main__":
	unittest.main()
