#!/usr/bin/env python3

import random
import unittest
from copy import copy

import networkit as nk

class TestMatchingAlgorithms(unittest.TestCase):

	def setUp(self):
		self.g = nk.readGraph("input/PGPgiantcompo.graph", nk.Format.METIS)
		self.gw = copy(self.g)
		nk.graphtools.randomizeWeights(self.gw)

	def hasUnmatchedNeighbors(self, g, m):
		for e in g.iterEdges():
			if not m.isMatched(e[0]) and not m.isMatched(e[1]):
				return True
		return False

	def testPathGrowingMatcher(self):
		def runAlgo(g):
			pgm = nk.matching.PathGrowingMatcher(self.g)
			pgm.run()
			m = pgm.getMatching()

		runAlgo(self.g)
		runAlgo(self.gw)

	def testSuitorMatcher(self):

		def doTest(g):
			m1 = nk.matching.SuitorMatcher(g, False).run().getMatching()
			nk.graphtools.sortEdgesByWeight(g, True)
			self.assertTrue(m1.isProper(g))
			self.assertFalse(self.hasUnmatchedNeighbors(g, m1))

			m2 = nk.matching.SuitorMatcher(g, True).run().getMatching()
			self.assertTrue(m2.isProper(g))
			self.assertFalse(self.hasUnmatchedNeighbors(g, m2))
			for u in g.iterNodes():
				self.assertEqual(m1.mate(u), m2.mate(u))

		doTest(self.g)
		doTest(self.gw)

if __name__ == "__main__":
	unittest.main()
