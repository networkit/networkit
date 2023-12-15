#!/usr/bin/env python3
import unittest

import networkit as nk


class TestSCD(unittest.TestCase):

	def setUp(self):
		self.G = nk.readGraph("input/PGPgiantcompo.graph", nk.Format.METIS)

	def testApproximatePageRank(self):
		epsilon = 1e-12
		result = dict(nk.scd.ApproximatePageRank(self.G, 0.1, epsilon).run(0))

		result2 = dict(nk.scd.ApproximatePageRank(self.G, 0.1, epsilon).run([0]))

		self.assertGreaterEqual(len(result), 1)
		self.assertEqual(len(result), len(result2))

		for u, score in result.items():
			self.assertTrue(self.G.hasNode(u))
			self.assertGreaterEqual(score, 0)
			self.assertTrue(u in result2)
			self.assertAlmostEqual(score, result2[u])

	def testSCD(self):
		nk.setSeed(42, False)
		seed = 20
		seeds = [seed]
		for name, algo in [
				("PageRankNibble", nk.scd.PageRankNibble(self.G, 0.1, 1e-12)),
				("GCE L", nk.scd.GCE(self.G, "L")),
				("GCE M", nk.scd.GCE(self.G, "M")),
				("LFM", nk.scd.LFMLocal(self.G, 0.8)),
				("TwoPhaseL", nk.scd.TwoPhaseL(self.G)),
				("TCE", nk.scd.TCE(self.G)),
				("LTE", nk.scd.LocalTightnessExpansion(self.G)),
				("LocalT", nk.scd.LocalT(self.G)),
				("Clique", nk.scd.CliqueDetect(self.G))]:
			result = algo.run(seeds)[seed]

			self.assertGreaterEqual(len(result), 1, "{} has empty community".format(name))

			cond = nk.scd.SetConductance(self.G, result).run().getConductance()

			self.assertLessEqual(cond, 0.5, "{} has too large conductance".format(name))

	def testSCDGroundTruth(self):
		G = nk.Graph(5, directed=True)
		G.addEdge(0,1)
		G.addEdge(2,1)
		G.addEdge(1,3)
		G.addEdge(0,4)
		G.addEdge(3,2)

		ground = nk.Cover(5)
		ground.allToSingletons()
		ground.addToSubset(1,2)
		ground.addToSubset(2,3)

		gt = nk.scd.SCDGroundTruthComparison(G, ground, {1:[1,2]}, True)
		gt.run()
		self.assertDictEqual(gt.getIndividualJaccard(), {1: 0.5})
		self.assertDictEqual(gt.getIndividualPrecision(), {1: 0.5})
		self.assertDictEqual(gt.getIndividualRecall(), {1: 1.0})
		self.assertDictEqual(gt.getIndividualF1(), {1: 0.6666666666666666})
		self.assertEqual(gt.getAveragePrecision(), 0.5)
		self.assertEqual(gt.getAverageRecall(), 1.0)
		self.assertAlmostEqual(gt.getAverageF1(), 0.667, 3)
		self.assertEqual(gt.getAverageJaccard(), 0.5)

	def testExpandOneCommunity(self):
		nk.setSeed(42, False)
		ground = nk.Cover(self.G.upperNodeIdBound())
		ground.setUpperBound(1)
		for i in range(10):
			ground.addToSubset(0, i)
		
		bfs = nk.scd.RandomBFS(self.G, ground)
		community = bfs.expandOneCommunity([10,20,30,40,50])	

		self.assertEqual(len(community), 5)
		self.assertRaises(TypeError, bfs.expandOneCommunity(5))

if __name__ == "__main__":
	unittest.main()
