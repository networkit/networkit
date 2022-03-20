#!/usr/bin/env python3
import unittest

import networkit as nk


class Test_SCD(unittest.TestCase):

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

if __name__ == "__main__":
	unittest.main()
