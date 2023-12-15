#!/usr/bin/env python3
import unittest
import os
import networkit as nk

class TestLinkprediction(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops
	
	def testAdamicAdarIndex(self):		
		AI = nk.linkprediction.AdamicAdarIndex()
		AI.setGraph(self.L)
		self.assertEqual(AI.run(0,7), 0.0)		

	def testAdjustedRandIndex(self):		
		ARI = nk.linkprediction.AdjustedRandIndex()
		ARI.setGraph(self.L)
		self.assertAlmostEqual(ARI.run(0,7), -0.44, delta=0.1)	

	def testAlgebraicDistanceIndex(self):		
		ADI = nk.linkprediction.AlgebraicDistanceIndex(self.L, 5, 10)
		ADI.setGraph(self.L)
		ADI.preprocess()
		self.assertAlmostEqual(ADI.run(0,7), 0.31, delta=0.25)	
	
	def testCommonNeighborsIndex(self):
		CNI = nk.linkprediction.CommonNeighborsIndex()
		CNI.setGraph(self.L)
		self.assertEqual(CNI.run(0,7), 0.0) 
		self.assertEqual(len(CNI.runAll()), 57)
		res = CNI.runOn([(0,1), (2,3), (3,7)])
		expectedRes = [((0, 1), 0.0), ((2, 3), 2.0), ((3, 7), 0.0)]
		for i in range(3):
			self.assertAlmostEqual(res[i][1], expectedRes[i][1], delta = 0.1)
	
	def testJaccardIndex(self):		
		JI = nk.linkprediction.JaccardIndex()
		JI.setGraph(self.L)
		self.assertEqual(JI.run(0,4), 0.25)	
	
	def testKatzIndex(self):
		KI = nk.linkprediction.KatzIndex()
		KI.setGraph(self.L)
		self.assertAlmostEqual(KI.run(0,7), 6.3125e-10, 2) 
		self.assertEqual(len(KI.runAll()), 57)
		res = KI.runOn([(0,1), (2,3), (3,7)])
		expectedRes = [((0, 1), 0.005), ((2, 3), 0.005), ((3, 7), 1.262e-07)]
		for i in range(3):
			self.assertAlmostEqual(res[i][1], expectedRes[i][1], delta = 0.1)

	def testLinkThresholder(self):
		LTH = nk.linkprediction.LinkThresholder()
		pred = [((0,1), 0.25), ((2,3), 2.5), ((4,7), 3.0)]		
		self.assertListEqual(LTH.byScore(pred, 0.3), [(2,3), (4,7)])
		self.assertListEqual(LTH.byCount(pred, 2), [(2,3), (4,7)])			

	def testMissingLinksFinder(self):
		MLF = nk.linkprediction.MissingLinksFinder(self.L)
		self.assertListEqual(MLF.findAtDistance(3),[(0, 6), (1, 5), (1, 7), (2, 5), (2, 7), (3, 5), (3, 7), (4, 8)])
		self.assertListEqual(MLF.findFromNode(0, 2), [(0, 2), (0, 3), (0, 4)])

	def testNeighborhoodDistanceIndex(self):		
		NDI = nk.linkprediction.NeighborhoodDistanceIndex()
		NDI.setGraph(self.L)
		self.assertEqual(NDI.run(0,7), 0.0)	

	def testNeighborhoodUtility(self):
		NHU = nk.linkprediction.NeighborhoodUtility()
		self.assertListEqual(NHU.getNeighborsUnion(self.L, 0, 5), [1,6,7])
		self.assertListEqual(NHU.getCommonNeighbors(self.L, 0, 3), [1])

	def testNeighborsMeasureIndex(self):		
		NMI = nk.linkprediction.NeighborsMeasureIndex()
		NMI.setGraph(self.L)
		self.assertEqual(NMI.run(0,7), 0.0)	

	def testPredictionSorter(self):
		PS = nk.linkprediction.PredictionsSorter()
		pred = [((0,1), 0.25), ((2,3), 12.5), ((4,7), 3.0)]		
		PS.sortByScore(pred)
		self.assertListEqual(pred, [((2,3), 12.5), ((4,7), 3.0), ((0,1), 0.25)])
		PS.sortByNodePair(pred)
		self.assertListEqual(pred, [((0,1), 0.25), ((2,3), 12.5), ((4,7), 3.0)])
		
	def testPreferentialAttachmentIndex(self):		
		PAI = nk.linkprediction.PreferentialAttachmentIndex()
		PAI.setGraph(self.L)
		self.assertEqual(PAI.run(0,7), 3.0)

	def testPrecisionRecallMetric(self):
		PRM = nk.linkprediction.PrecisionRecallMetric(self.L)
		PRM.setTestGraph(self.L)
		pred = [((0,1),5), ((2,3), 2.5), ((4,7), 3)]
		self.assertListEqual(PRM.getCurve(pred)[0], [0.0, 0.5, 1.0])
		self.assertEqual(PRM.getCurve(pred)[1][0], 1.0)
		self.assertEqual(PRM.getCurve(pred)[1][1], 0.5)
		self.assertAlmostEqual(PRM.getCurve(pred)[1][2], 0.667, 3)
		self.assertAlmostEqual(PRM.getAreaUnderCurve(), 0.667 , 3)

	def testRandomLinkSampler(self):
		RLS = nk.linkprediction.RandomLinkSampler(self.L, 50)
		resGraph = (RLS.byCount(self.L, 3))
		self.assertEqual(resGraph.numberOfEdges(), 3)

	def testROCMetric(self):
		RM = nk.linkprediction.ROCMetric(self.L)
		RM.setTestGraph(self.L)
		pred = [((0,1),5), ((2,3), 2.5), ((4,7), 3)]
		self.assertTupleEqual(RM.getCurve(pred), ([0.0, 1.0], [0.5, 1.0]))	
		self.assertEqual(RM.getAreaUnderCurve(), 0.75)	

	def testResourceAllocationIndex(self):		
		I = nk.linkprediction.ResourceAllocationIndex()
		I.setGraph(self.L)
		self.assertEqual(I.run(0,4), 0.25)	

	def testSameCommunityIndex(self):		
		SCI = nk.linkprediction.SameCommunityIndex()
		SCI.setGraph(self.L)
		self.assertEqual(SCI.run(0,4), 1.0)
	
	def testTotalNeighborsIndex(self):		
		TNI = nk.linkprediction.TotalNeighborsIndex()
		TNI.setGraph(self.L)
		self.assertEqual(TNI.run(0,7), 4.0)	

	def testUDegreeIndex(self):		
		UI = nk.linkprediction.UDegreeIndex()
		UI.setGraph(self.L)
		self.assertEqual(UI.run(0,7), 1.0)

	def testVDegreeIndex(self):		
		VI = nk.linkprediction.VDegreeIndex()
		VI.setGraph(self.L)
		self.assertEqual(VI.run(0,7), 3.0)													

if __name__ == "__main__":
	unittest.main()

