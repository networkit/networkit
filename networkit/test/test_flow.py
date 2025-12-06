#!/usr/bin/env python3
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

	def testDinicConstructorPreconditions(self):
		# Undirected graph
		G_undirected = nk.Graph(2, weighted=True, directed=False)
		with self.assertRaises(RuntimeError) as ctx1:
			nk.flow.Dinic(G_undirected, 0, 1)
		self.assertEqual(str(ctx1.exception), "Dinic algorithm requires directed graph!")

		# Unweighted graph
		G_unweighted = nk.Graph(2, weighted=False, directed=True)
		with self.assertRaises(RuntimeError) as ctx2:
			nk.flow.Dinic(G_unweighted, 0, 1)
		self.assertEqual(str(ctx2.exception), "Dinic algorithm requires weighted graph!")

		# Same source and target
		G_same = nk.Graph(2, weighted=True, directed=True)
		with self.assertRaises(RuntimeError) as ctx3:
			nk.flow.Dinic(G_same, 0, 0)
		self.assertEqual(
			str(ctx3.exception),
			"Dinic algorithm requires `source` and `target` node to be different!",
		)

	def testDinicRunNotCalledThrows(self):
		G = nk.Graph(2, weighted=True, directed=True)
		algo = nk.flow.Dinic(G, 0, 1)
		with self.assertRaises(RuntimeError) as ctx:
			algo.getMaxFlow()
		self.assertEqual(str(ctx.exception), "Error, run must be called first")

	def testDinicNegativeWeightsThrows(self):
		G = nk.Graph(3, weighted=True, directed=True)
		G.addEdge(0, 1, 1.0)
		G.addEdge(1, 2, -1.0)
		algo = nk.flow.Dinic(G, 0, 2)
		with self.assertRaises(RuntimeError):
			algo.run()

	def testDinicThreeDisjointPathsWithParallelEdges(self):
		"""
		Combines:
		- multiple edge-disjoint s–t paths
		- parallel edges + indexEdges()
		Expected max flow: 3 (three unit-capacity branches)
		"""
		G = nk.Graph(5, weighted=True, directed=True)
		# branch via 1
		G.addEdge(0, 1, 1.0)
		G.addEdge(1, 0, 1.0)  # back edge
		G.addEdge(1, 4, 1.0)
		# branch via 2
		G.addEdge(0, 2, 1.0)
		G.addEdge(2, 0, 1.0)  # back edge
		G.addEdge(2, 4, 1.0)
		# branch via 3
		G.addEdge(0, 3, 1.0)
		G.addEdge(3, 0, 1.0)  # back edge
		G.addEdge(3, 4, 1.0)
		G.indexEdges()

		algo = nk.flow.Dinic(G, 0, 4)
		algo.run()
		self.assertAlmostEqual(algo.getMaxFlow(), 3.0, places=12)

	def testDinicDisconnectedGraphs(self):
		"""
		Source and sink in different connected components → max flow must be zero.
		"""
		G = nk.Graph(7, weighted=True, directed=True)
		# component 1
		G.addEdge(0, 1, 10.0)
		G.addEdge(1, 2, 5.0)
		G.addEdge(2, 3, 7.0)
		# component 2
		G.addEdge(4, 5, 11.0)
		G.addEdge(5, 6, 10.0)

		algo = nk.flow.Dinic(G, 0, 5)
		algo.run()
		self.assertAlmostEqual(algo.getMaxFlow(), 0.0, places=12)

	def testDinicNumericalStabilityDecimalSplits(self):
		"""
		Parallel s–t paths with small fractional capacities
		plus a sub-epsilon direct edge that should be ignored
		by the tolerance gating.
		Total incoming capacity at sink should be 1.0.
		"""
		G = nk.Graph(7, weighted=True, directed=True)
		G.addEdge(0, 1, 1.0)
		G.addEdge(1, 2, 0.1)
		G.addEdge(2, 6, 0.1)
		G.addEdge(1, 3, 0.2)
		G.addEdge(3, 6, 0.2)
		G.addEdge(1, 4, 0.3)
		G.addEdge(4, 6, 0.3)
		G.addEdge(1, 5, 0.4)
		G.addEdge(5, 6, 0.4)

		# sub-epsilon direct edge that must effectively be ignored
		G.addEdge(0, 6, 1e-18)

		algo = nk.flow.Dinic(G, 0, 6)
		algo.run()
		self.assertAlmostEqual(algo.getMaxFlow(), 1.0, places=12)

if __name__ == "__main__":
	unittest.main()
