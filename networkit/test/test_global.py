#!/usr/bin/env python3
import numpy as np
import unittest

import networkit as nk

class TestGlobal(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops
		self.LL = nk.readGraph("input/looptest2.gml", nk.Format.GML) #with self-loops sprinkled in
	
	def testGlobalsClusteringCoefficient(self):
		CL = nk.globals.ClusteringCoefficient()
		CL.exactGlobal(self.L)
		CL.exactGlobal(self.LL)
		CL.approxGlobal(self.L, 5)
		CL.approxGlobal(self.LL, 5)
		CL.approxAvgLocal(self.L, 5)
		CL.approxAvgLocal(self.LL, 5)
		CL.sequentialAvgLocal(self.L)
		CL.sequentialAvgLocal(self.LL)

	def testStatsGiniCoefficient(self):
		# Perfect equality -> 0
		self.assertAlmostEqual(nk.stats.gini([1, 1, 1, 1]), 0.0, places=12)
		# permutation invariance
		self.assertAlmostEqual(
			nk.stats.gini([3, 1, 2]),
			nk.stats.gini([1, 2, 3]),
			places=12
		)
		# scale invariance
		self.assertAlmostEqual(
			nk.stats.gini([1, 2, 3]),
			nk.stats.gini([10, 20, 30]),
			places=12
		)
		# expected 0.75
		self.assertAlmostEqual(nk.stats.gini([0, 0, 0, 4]), 0.75, places=12)
		# numpy array as input
		arr = np.array([0, 0, 0, 4], dtype=np.float64)
		self.assertAlmostEqual(nk.stats.gini(arr), 0.75, places=12)

if __name__ == "__main__":
	unittest.main()
