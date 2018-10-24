#! /usr/bin/python

import unittest
import os

from networkit import *

class Test_SelfLoops(unittest.TestCase):

	def checkCovers(self, c1, c2):
		if not c1.numberOfElements() == c2.numberOfElements(): return False
		if not c1.numberOfSubsets() == c2. numberOfSubsets(): return False
		for i in range(0,c1.numberOfElements()):
			if not c1.subsetsOf(i) == c2.subsetsOf(i): return False
		return True

	def setUp(self):
		# toggle the comment/uncomment to test on small or large test cases
		#self.L = readGraph("PGPgiantcompo.graph", Format.METIS) #without self-loops
		#self.LL = readGraph("PGPConnectedCompoLoops.gml", Format.GML) #with self-loops sprinkled in
		self.L = readGraph("input/looptest1.gml", Format.GML) #without self-loops
		self.LL = readGraph("input/looptest2.gml", Format.GML) #with self-loops sprinkled in

	def test_centrality_Betweenness(self):
		CL = centrality.Betweenness(self.L)
		CL.run()
		CLL = centrality.Betweenness(self.LL)
		CLL.run()
		self.assertEqual(CL.ranking(), CLL.ranking())

	def test_centrality_ApproxBetweenness(self):
		CL = centrality.ApproxBetweenness(self.L, epsilon=0.01, delta=0.1)
		CL.run()
		CLL = centrality.ApproxBetweenness(self.LL, epsilon=0.01, delta=0.1)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))
		for i in range(len(CL.ranking())):
			self.assertAlmostEqual(CL.ranking()[i][1], CLL.ranking()[i][1], delta=0.2*CL.ranking()[i][1])


	def test_centrality_Closeness(self):
		CL = centrality.Closeness(self.L, True, centrality.ClosenessVariant.Generalized)
		CL.run()
		CLL = centrality.Closeness(self.LL, True, centrality.ClosenessVariant.Generalized)
		CLL.run()
		self.assertEqual(CL.ranking(), CLL.ranking())


	def test_centrality_CoreDecomposition(self):
		CL = centrality.CoreDecomposition(self.L)
		CL.run()
		try:
			CLL = centrality.CoreDecomposition(self.LL)
		except RuntimeError:
			import copy
			tmp = copy.deepcopy(self.LL)
			tmp.removeSelfLoops()
			CLL = centrality.CoreDecomposition(tmp)
			CLL.run()
			self.assertTrue(self.checkCovers(CL.getCover(),CLL.getCover()))


	def test_centrality_EigenvectorCentrality(self):
		CL = centrality.EigenvectorCentrality(self.L)
		CL.run()
		CLL = centrality.EigenvectorCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def test_centrality_KPathCentrality(self):
		CL = centrality.KPathCentrality(self.L)
		CL.run()
		CLL = centrality.KPathCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def test_centrality_KatzCentrality(self):
		CL = centrality.KatzCentrality(self.L)
		CL.run()
		CLL = centrality.KatzCentrality(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def test_centrality_PageRank(self):
		CL = centrality.PageRank(self.L)
		CL.run()
		CLL = centrality.PageRank(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def test_centrality_rankPerNode(self):
		CL = centrality.PageRank(self.L)
		CL.run()
		CLL = centrality.PageRank(self.LL)
		CLL.run()
		#test if list of pairs and list of ranks have the same length
		self.assertEqual(len(CL.ranking()),len(centrality.rankPerNode(CL.ranking())))
		self.assertEqual(len(CLL.ranking()),len(centrality.rankPerNode(CLL.ranking())))


	def test_centrality_SciPyPageRank(self):
		CL = centrality.SciPyPageRank(self.L)
		CL.run()
		CLL = centrality.SciPyPageRank(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))


	def test_centrality_SciPyEVZ(self):
		CL = centrality.SciPyEVZ(self.L)
		CL.run()
		CLL = centrality.SciPyEVZ(self.LL)
		CLL.run()
		#test if lists have the same length
		self.assertEqual(len(CL.ranking()),len(CLL.ranking()))

	def test_centrality_relativeRankErrors(self):
		CL = centrality.Betweenness(self.L)
		CL.run()
		CLL = centrality.Betweenness(self.LL)
		CLL.run()
		self.assertEqual(len(CL.ranking()), len(centrality.relativeRankErrors(CL.ranking(),CLL.ranking())))


	def test_community_PLM(self):
		PLML = community.PLM(self.L)
		PLMLL = community.PLM(self.LL)
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
		self.assertEqual(set(self.L.nodes()), set(reconstructedSet))
		reconstructedSet = []
		for i in PLMLLP.getSubsetIds():
			for j in PLMLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.LL.nodes()), set(reconstructedSet))

	def test_community_PLP(self):
		PLPL = community.PLP(self.L)
		PLPLL = community.PLP(self.LL)
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
		self.assertEqual(set(self.L.nodes()), set(reconstructedSet))
		reconstructedSet = []
		for i in PLPLLP.getSubsetIds():
			for j in PLPLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.LL.nodes()), set(reconstructedSet))


	def test_community_CutClustering(self):
		CL = community.CutClustering(self.L, 0.2)
		CLL = community.CutClustering(self.LL, 0.2)
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
		self.assertEqual(set(self.L.nodes()), set(reconstructedSet))
		reconstructedSet = []
		for i in CLLP.getSubsetIds():
			for j in CLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.LL.nodes()), set(reconstructedSet))


	def test_community_GraphClusteringTools(self):
		PLMLL = community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		GCT = community.GraphClusteringTools()
		self.assertIsInstance(GCT.equalClustering(PLMLLP,PLPLLP, self.LL), bool)
		self.assertIsInstance(GCT.getImbalance(PLMLLP), float)
		self.assertIsInstance(GCT.isOneClustering(self.LL, PLPLLP), bool)
		self.assertIsInstance(GCT.isProperClustering(self.LL, PLMLLP), bool)
		self.assertIsInstance(GCT.isSingletonClustering(self.LL, PLPLLP), bool)


	def test_community_GraphStructuralRandMeasure(self):
		PLMLL = community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		GSRM = community.GraphStructuralRandMeasure()
		self.assertAlmostEqual(GSRM.getDissimilarity(self.LL, PLMLLP, PLPLLP),0.5, delta=0.5 )


	def test_community_Hubdominance(self):
		PLMLL = community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		HD = community.HubDominance()
		self.assertIsInstance(HD.getQuality(PLMLLP, self.LL),float )


	def test_community_JaccardMeasure(self):
		PLMLL = community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		JM = community.JaccardMeasure()
		self.assertIsInstance(JM.getDissimilarity(self.LL, PLMLLP, PLPLLP),float)


	def test_community_LPDegreeOrdered(self):
		CL = community.LPDegreeOrdered(self.L)
		CLL = community.LPDegreeOrdered(self.LL)
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
		self.assertEqual(set(self.L.nodes()), set(reconstructedSet))
		reconstructedSet = []
		for i in CLLP.getSubsetIds():
			for j in CLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.LL.nodes()), set(reconstructedSet))

	def test_community_Modularity(self):
		PLPLL = community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		Mod = community.Modularity()
		self.assertAlmostEqual(Mod.getQuality(PLPLLP, self.LL),0.25, delta=0.75)


	def test_community_NMIDistance(self):
		PLMLL = community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		NMI = community.NMIDistance()
		self.assertIsInstance(NMI.getDissimilarity(self.LL, PLMLLP, PLPLLP),float)


	def test_community_NodeStructuralRandMeasure(self):
		PLMLL = community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		NSRM = community.NodeStructuralRandMeasure()
		self.assertAlmostEqual(NSRM.getDissimilarity(self.LL, PLMLLP, PLPLLP),0.5, delta=0.5 )


	def test_community_communityGraph(self):
		PLMLL = community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		CG = community.communityGraph(self.LL, PLMLLP)
		self.assertIsInstance(len(CG.nodes()), int)


	def test_community_evaluateCommunityDetection(self):
		PLMLL = community.PLM(self.LL)
		community.evalCommunityDetection(PLMLL, self.LL)


	def test_community_kCoreCommunityDetection(self):
		with self.assertRaises(RuntimeError) as cm:
			kCCD = community.kCoreCommunityDetection(self.LL, 1, inspect=False)


	def test_flow_EdmondsKarp(self):
		self.L.indexEdges()
		self.LL.indexEdges()
		r1 = self.L.randomNode()
		r2 = self.L.randomNode()
		while r1 is r2:
			r2 = self.L.randomNode()
		EKL = flow.EdmondsKarp(self.L, r1, r2)
		EKLL = flow.EdmondsKarp(self.LL, r1, r2)
		EKL.run()
		EKLL.run()


	def test_clique_MaxClique(self):
		clique.MaxClique(self.LL).run()


	def test_globals_ClusteringCoefficient(self):
		CL = globals.ClusteringCoefficient()
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


	def test_components_ConnectedComponents(self):
		CC = components.ConnectedComponents(self.LL)
		CC.run()
		CC.componentOfNode(1)
		CC.getComponentSizes()
		CC.getPartition()
		CC.numberOfComponents()


	def test_distance_Diameter(self):
		D = distance.Diameter(self.LL, distance.DiameterAlgo.EstimatedRange, error = 0.1)
		D.run()
		D = distance.Diameter(self.LL, distance.DiameterAlgo.EstimatedSamples, nSamples = 5)
		D.run()
		D = distance.Diameter(self.LL, distance.DiameterAlgo.Exact)
		D.run()


	def test_distance_Eccentricity(self):
		E = distance.Eccentricity()
		E.getValue(self.LL, 0)


	def test_distance_EffectiveDiameter(self):
		algo = distance.EffectiveDiameter(self.L)
		algo.run()
		algo = distance.EffectiveDiameter(self.LL)
		algo.run()


	def test_distance_ApproxEffectiveDiameter(self):
		algo = distance.EffectiveDiameterApproximation(self.L)
		algo.run()
		algo = distance.EffectiveDiameterApproximation(self.LL)
		algo.run()


	def test_distance_ApproxHopPlot(self):
		algo = distance.HopPlotApproximation(self.L)
		algo.run()
		algo = distance.HopPlotApproximation(self.LL)
		algo.run()


	def test_distance_NeighborhoodFunction(self):
		algo = distance.NeighborhoodFunction(self.L)
		algo.run()
		algo = distance.NeighborhoodFunction(self.LL)
		algo.run()


	def test_distance_ApproxNeighborhoodFunction(self):
		algo = distance.NeighborhoodFunctionApproximation(self.L)
		algo.run()
		algo = distance.NeighborhoodFunctionApproximation(self.LL)
		algo.run()


if __name__ == "__main__":
	unittest.main()
