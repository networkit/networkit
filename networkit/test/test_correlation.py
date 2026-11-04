#!/usr/bin/env python3
import unittest
import os
import networkit as nk

class TestCorrelation(unittest.TestCase):

	def testAssortativity(self):
		G=nk.Graph(4)
		G.addEdge(0,1)
		G.addEdge(2,3)
		assor = nk.correlation.Assortativity(G, [0.2, 0.4, 1.2, 2.2])
		assor.run()
		self.assertEqual(assor.getCoefficient(), 1.0)

	def testAssortativityPartition(self):
		G=nk.Graph(4)
		G.addEdge(0,1)
		G.addEdge(2,3)
		G.addEdge(1,2)
		P = nk.Partition(4)
		P.setUpperBound(2)
		P.addToSubset(0,0)
		P.addToSubset(0,1)
		P.addToSubset(1,2)
		P.addToSubset(1,3)
		assor = nk.correlation.Assortativity(G, P)
		assor.run()
		self.assertAlmostEqual(assor.getCoefficient(), 1 / 3)

	def testAssortativityPartitionSelfLoops(self):
		for directed in (False, True):
			G = nk.Graph(2, directed=directed)
			G.addEdge(0, 0)
			G.addEdge(0, 1)
			G.addEdge(1, 1)
			P = nk.Partition(2)
			P.setUpperBound(2)
			P.addToSubset(0, 0)
			P.addToSubset(1, 1)
			assor = nk.correlation.Assortativity(G, P).run()
			self.assertAlmostEqual(assor.getCoefficient(), 2 / 5 if directed else 0)

if __name__ == "__main__":
	unittest.main()
