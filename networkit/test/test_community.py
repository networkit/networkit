#!/usr/bin/env python3
import unittest
import os
import networkit as nk 


class TestCommunity(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops
		self.LL = nk.readGraph("input/looptest2.gml", nk.Format.GML) #with self-loops	

	def testCutClustering(self):
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

	def testGraph(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		CG = nk.community.communityGraph(self.LL, PLMLLP)

	def testEvaluateCommunityDetection(self):
		PLMLL = nk.community.PLM(self.LL)
		nk.community.evalCommunityDetection(PLMLL, self.LL)


	def testKCoreCommunityDetection(self):
		with self.assertRaises(RuntimeError) as cm:
			kCCD = nk.community.kCoreCommunityDetection(self.LL, 1, inspect=False)	
	
	def testGraphClusteringTools(self):
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

	def testGraphStructuralRandMeasure(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		GSRM = nk.community.GraphStructuralRandMeasure()
		self.assertAlmostEqual(GSRM.getDissimilarity(self.LL, PLMLLP, PLPLLP),0.5, delta=0.5 )

	def testHubdominance(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		HD = nk.community.HubDominance()
		self.assertIsInstance(HD.getQuality(PLMLLP, self.LL),float )


	def testJaccardMeasure(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		JM = nk.community.JaccardMeasure()
		self.assertIsInstance(JM.getDissimilarity(self.LL, PLMLLP, PLPLLP),float)

	def testLPDegreeOrdered(self):
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

	def testNMIDistance(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		NMI = nk.community.NMIDistance()
		self.assertIsInstance(NMI.getDissimilarity(self.LL, PLMLLP, PLPLLP),float)

	def testNodeStructuralRandMeasure(self):
		PLMLL = nk.community.PLM(self.LL)
		PLMLL.run()
		PLMLLP = PLMLL.getPartition()
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		NSRM = nk.community.NodeStructuralRandMeasure()
		self.assertAlmostEqual(NSRM.getDissimilarity(self.LL, PLMLLP, PLPLLP),0.5, delta=0.5 )		

	def testModularity(self):
		PLPLL = nk.community.PLP(self.LL)
		PLPLL.run()
		PLPLLP = PLPLL.getPartition()
		Mod = nk.community.Modularity()
		self.assertAlmostEqual(Mod.getQuality(PLPLLP, self.LL),0.25, delta=0.75)	

	def testPLM(self):
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

	def testPL(self):
		PLL = nk.community.ParallelLeiden(self.L)
		PLLL = nk.community.ParallelLeiden(self.LL)
		PLL.run()
		PLLL.run()
		PLLP = PLL.getPartition()
		PLLLP = PLLL.getPartition()
		self.assertIsNot(PLLP.getSubsetIds(), None)
		self.assertIsNot(PLLLP.getSubsetIds(), None)
		# test if partitions add up to original set
		reconstructedSet = []
		for i in PLLP.getSubsetIds():
			for j in PLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.L.iterNodes()), set(reconstructedSet))
		reconstructedSet = []
		for i in PLLLP.getSubsetIds():
			for j in PLLLP.getMembers(i):
				reconstructedSet.append(j)
		self.assertEqual(set(self.LL.iterNodes()), set(reconstructedSet))

	def testPLP(self):
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
	
	def testCutClustering(self):
		jazz = nk.readGraph("input/jazz.graph",nk.Format.METIS)

		cc = nk.community.CutClustering(jazz,0.5076142131979697)

		comms = nk.community.detectCommunities(jazz, cc)
		self.assertEqual(193, comms.numberOfSubsets())

		hierarchy = nk.community.CutClustering.getClusterHierarchy(jazz)
		self.assertEqual(3, len(hierarchy))
	
	def testLFM(self):
		G = nk.readGraph("input/PGPgiantcompo.graph",nk.Format.METIS)
		local = nk.scd.LFMLocal(G)
		lfm = nk.community.LFM(G, local)
		lfm.run()
		cover = lfm.getCover()
		for u in G.iterNodes():
			self.assertGreaterEqual(len(cover.subsetsOf(u)), 1)
	
	def testPLM(self):
		G = nk.readGraph("input/PGPgiantcompo.graph",nk.Format.METIS)
		comms = nk.community.detectCommunities(G, nk.community.PLM(G))

		# each node must be assigned
		for v in G.iterNodes():
			self.assertTrue(comms.contains(v))			
	

if __name__ == "__main__":
	unittest.main()