#!/usr/bin/env python3
import numpy as np
import os
import random
import unittest

import networkit as nk

class TestFlow(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops
		self.LL = nk.readGraph("input/looptest2.gml", nk.Format.GML) #with self-loops sprinkled in

	def testEdmondsKarp(self):
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

if __name__ == "__main__":
	unittest.main()
