#!/usr/bin/env python3
import numpy as np
import os
import random
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
		CL.avgLocal(self.L)
		with self.assertRaises(RuntimeError) as cm:
			CL.avgLocal(self.LL)
		CL.sequentialAvgLocal(self.L)
		CL.sequentialAvgLocal(self.LL)

if __name__ == "__main__":
	unittest.main()
