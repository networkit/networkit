#!/usr/bin/env python3
import numpy as np
import os
import random
import unittest

import networkit as nk

class TestFlow(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops

	def testEdmondsKarp(self):
		G = nk.Graph(5)
		G.addEdge(0,1)
		G.addEdge(1,2)
		G.addEdge(1,3)
		G.addEdge(3,4)
		G.addEdge(2,3)
		G.indexEdges()
		EKL = nk.flow.EdmondsKarp(G, 0, 1)
		EKL.run()
		self.assertEqual(EKL.getFlow(0, 1), 1.0)
		self.assertEqual(EKL.getFlow(1), 0.0)
		self.assertEqual(EKL.getMaxFlow(), 1.0)
		self.assertListEqual(EKL.getSourceSet(), [0])
		self.assertListEqual(EKL.getFlowVector(), [-1.0, 0.0, 0.0, 0.0, 0.0])

if __name__ == "__main__":
	unittest.main()
