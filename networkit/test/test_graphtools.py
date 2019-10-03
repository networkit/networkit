#!/usr/bin/env python3
import unittest
import random
import networkit as nk

class TestGraphTools(unittest.TestCase):
	def getSmallGraph(self, weighted=False, directed=False):
		G = nk.Graph(4, weighted, directed)
		G.addEdge(0, 1, 1.0)
		G.addEdge(0, 2, 2.0)
		G.addEdge(3, 1, 4.0)
		G.addEdge(3, 2, 5.0)
		G.addEdge(1, 2, 3.0)

		return G

	def generateRandomWeights(self, G):
		if not G.isWeighted():
			G = nk.graph.GraphTools.toWeighted(G)
		G.forEdges(lambda u, v, w, eid: G.setWeight(u, v, random.random()))

		return G

	def testCopyNodes(self):
		def checkNodes(G, GCopy):
			self.assertEqual(G.isDirected(), GCopy.isDirected())
			self.assertEqual(G.isWeighted(), GCopy.isWeighted())
			self.assertEqual(GCopy.numberOfEdges(), 0)
			for u in range(G.upperNodeIdBound()):
				self.assertEqual(G.hasNode(u), GCopy.hasNode(u))

		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)
				GCopy = nk.graph.GraphTools.copyNodes(G)
				checkNodes(G, GCopy)

				for _ in range(1, G.numberOfNodes()):
					G.removeNode(G.randomNode())
					GCopy = nk.graph.GraphTools.copyNodes(G)
					checkNodes(G, GCopy)

	def testSubgraphFromNodesUndirected(self):
		G = self.getSmallGraph(True, False)

		nodes = set([0])
		res = nk.graph.GraphTools.subgraphFromNodes(G, nodes)
		self.assertTrue(res.isWeighted())
		self.assertFalse(res.isDirected())
		self.assertEqual(res.numberOfNodes(), 1)
		self.assertEqual(res.numberOfEdges(), 0)

		res = nk.graph.GraphTools.subgraphFromNodes(G, nodes, True)

		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 2) # 0-1, 0-2, NOT 1-2

		self.assertEqual(G.weight(0, 1), 1.0)
		self.assertEqual(G.weight(0, 2), 2.0)

		nodes = set([0, 1])
		res = nk.graph.GraphTools.subgraphFromNodes(G, nodes)
		self.assertEqual(res.numberOfNodes(), 2)
		self.assertEqual(res.numberOfEdges(), 1) # 0 - 1

		res = nk.graph.GraphTools.subgraphFromNodes(G, nodes, True)
		self.assertEqual(res.numberOfNodes(), 4)
		self.assertEqual(res.numberOfEdges(), 4) # 0-1, 0-2, 1-2, 1-3

	def testSubgraphFromNodesDirected(self):
		G = self.getSmallGraph(True, True)

		nodes = set([0])
		res = nk.graph.GraphTools.subgraphFromNodes(G, nodes)

		self.assertTrue(res.isWeighted())
		self.assertTrue(res.isDirected())

		self.assertEqual(res.numberOfNodes(), 1)
		self.assertEqual(res.numberOfEdges(), 0)

		nodes = set([0])
		res = nk.graph.GraphTools.subgraphFromNodes(G, nodes, True)
		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 2) # 0->1, 0->2, NOT 1->2

		nodes = set([0, 1])
		res = nk.graph.GraphTools.subgraphFromNodes(G, nodes)
		self.assertEqual(res.numberOfNodes(), 2)
		self.assertEqual(res.numberOfEdges(), 1) # 0 -> 1

		nodes = set([0, 1])
		res = nk.graph.GraphTools.subgraphFromNodes(G, nodes, True)
		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 3) # 0->1, 0->2, 1->2

		nodes = set([0, 1])
		res = nk.graph.GraphTools.subgraphFromNodes(G, nodes, True, True)
		self.assertEqual(res.numberOfNodes(), 4)
		self.assertEqual(res.numberOfEdges(), 4) # 0->1, 0->2, 1->2, 3->1

	def testGraphTranspose(self):
		for seed in range(1, 4):
			nk.setSeed(seed, True)
			random.seed(seed)
			G = nk.generators.ErdosRenyiGenerator(100, 0.2, True).generate()

			for _ in range(20):
				u = G.randomNode()
				if not G.hasEdge(u, u):
					G.addEdge(u, u)

			# Delete a few nodes
			for _ in range(10):
				G.removeNode(G.randomNode())
			self.assertGreater(G.numberOfSelfLoops(), 0)

			# Assign random weights
			GWeighted = self.generateRandomWeights(G)

			GWeighted.indexEdges()
			GTrans = nk.graph.GraphTools.transpose(GWeighted)

			def checkGWeightedEdges(u, v, w, eid):
				self.assertEqual(GWeighted.edgeId(u, v), GTrans.edgeId(v, u))
				self.assertEqual(GWeighted.weight(u, v), GTrans.weight(v, u))
			GWeighted.forEdges(checkGWeightedEdges)

			def checkGTransEdges(v, u, w, eid):
				self.assertEqual(GWeighted.edgeId(u, v), GTrans.edgeId(v, u))
				self.assertEqual(GWeighted.weight(u, v), GTrans.weight(v, u))
			GTrans.forEdges(checkGTransEdges)

			for u in range(GWeighted.upperNodeIdBound()):
				self.assertEqual(GWeighted.hasNode(u), GTrans.hasNode(u))

			self.assertEqual(GWeighted.numberOfNodes(), GTrans.numberOfNodes())
			self.assertEqual(GWeighted.upperNodeIdBound(), GTrans.upperNodeIdBound())
			self.assertEqual(GWeighted.numberOfEdges(), GTrans.numberOfEdges())
			self.assertEqual(GWeighted.upperEdgeIdBound(), GTrans.upperEdgeIdBound())
			self.assertEqual(GWeighted.numberOfSelfLoops(), GTrans.numberOfSelfLoops())

	def testToUndirected(self):
		n = 200
		p = 0.2

		def testGraphs(G, G1):
			self.assertEqual(G.numberOfNodes(), G1.numberOfNodes())
			self.assertEqual(G.upperNodeIdBound(), G1.upperNodeIdBound())
			self.assertEqual(G.numberOfEdges(), G1.numberOfEdges())
			self.assertEqual(G.upperEdgeIdBound(), G1.upperEdgeIdBound())
			self.assertEqual(G.isWeighted(), G1.isWeighted())
			self.assertNotEqual(G.isDirected(), G1.isDirected())
			self.assertEqual(G.hasEdgeIds(), G1.hasEdgeIds())

			def testEdges(u, v, w, eid):
				self.assertTrue(G1.hasEdge(u, v))
				self.assertEqual(G1.weight(u, v), w)
			G.forEdges(testEdges)

		for seed in range(1, 4):
			nk.setSeed(seed, False)
			random.seed(seed)
			G = nk.generators.ErdosRenyiGenerator(n, p, True).generate()
			for weighted in [True, False]:
				if weighted:
					G = self.generateRandomWeights(G)
				G1 = nk.graph.GraphTools.toUndirected(G)
				testGraphs(G, G1)

	def testToUnWeighted(self):
		n = 200
		p = 0.2

		def testGraphs(G, G1):
			self.assertEqual(G.numberOfNodes(), G1.numberOfNodes())
			self.assertEqual(G.upperNodeIdBound(), G1.upperNodeIdBound())
			self.assertEqual(G.numberOfEdges(), G1.numberOfEdges())
			self.assertNotEqual(G.isWeighted(), G1.isWeighted())
			self.assertEqual(G.isDirected(), G1.isDirected())
			self.assertEqual(G.hasEdgeIds(), G1.hasEdgeIds())

			def checkEdges(u, v, w, eid):
				self.assertTrue(G1.hasEdge(u, v))
				if G1.isWeighted():
					self.assertEqual(G1.weight(u, v), 1.0)
			G.forEdges(checkEdges)

		for seed in range(1, 4):
			nk.setSeed(seed, False)
			random.seed(seed)
			for directed in [True, False]:
				G = nk.generators.ErdosRenyiGenerator(n, p, directed).generate()

				G1 = nk.graph.GraphTools.toWeighted(G)
				testGraphs(G, G1)

				G = self.generateRandomWeights(G)

				G1 = nk.graph.GraphTools.toUnweighted(G)
				testGraphs(G, G1)

if __name__ == "__main__":
	unittest.main()
