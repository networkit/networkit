#!/usr/bin/env python3
import unittest
import random
from copy import copy
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
			G = nk.graphtools.toWeighted(G)
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
				GCopy = nk.graphtools.copyNodes(G)
				checkNodes(G, GCopy)

				for _ in range(1, G.numberOfNodes()):
					G.removeNode(G.randomNode())
					GCopy = nk.graphtools.copyNodes(G)
					checkNodes(G, GCopy)

	def testSubgraphFromNodesUndirected(self):
		G = self.getSmallGraph(True, False)

		nodes = set([0])
		res = nk.graphtools.subgraphFromNodes(G, nodes)
		self.assertTrue(res.isWeighted())
		self.assertFalse(res.isDirected())
		self.assertEqual(res.numberOfNodes(), 1)
		self.assertEqual(res.numberOfEdges(), 0)

		res = nk.graphtools.subgraphFromNodes(G, nodes, True)

		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 2) # 0-1, 0-2, NOT 1-2

		self.assertEqual(G.weight(0, 1), 1.0)
		self.assertEqual(G.weight(0, 2), 2.0)

		nodes = set([0, 1])
		res = nk.graphtools.subgraphFromNodes(G, nodes)
		self.assertEqual(res.numberOfNodes(), 2)
		self.assertEqual(res.numberOfEdges(), 1) # 0 - 1

		res = nk.graphtools.subgraphFromNodes(G, nodes, True)
		self.assertEqual(res.numberOfNodes(), 4)
		self.assertEqual(res.numberOfEdges(), 4) # 0-1, 0-2, 1-2, 1-3

	def testSubgraphFromNodesDirected(self):
		G = self.getSmallGraph(True, True)

		nodes = set([0])
		res = nk.graphtools.subgraphFromNodes(G, nodes)

		self.assertTrue(res.isWeighted())
		self.assertTrue(res.isDirected())

		self.assertEqual(res.numberOfNodes(), 1)
		self.assertEqual(res.numberOfEdges(), 0)

		nodes = set([0])
		res = nk.graphtools.subgraphFromNodes(G, nodes, True)
		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 2) # 0->1, 0->2, NOT 1->2

		nodes = set([0, 1])
		res = nk.graphtools.subgraphFromNodes(G, nodes)
		self.assertEqual(res.numberOfNodes(), 2)
		self.assertEqual(res.numberOfEdges(), 1) # 0 -> 1

		nodes = set([0, 1])
		res = nk.graphtools.subgraphFromNodes(G, nodes, True)
		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 3) # 0->1, 0->2, 1->2

		nodes = set([0, 1])
		res = nk.graphtools.subgraphFromNodes(G, nodes, True, True)
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
			GTrans = nk.graphtools.transpose(GWeighted)

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
				G1 = nk.graphtools.toUndirected(G)
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

				G1 = nk.graphtools.toWeighted(G)
				testGraphs(G, G1)

				G = self.generateRandomWeights(G)

				G1 = nk.graphtools.toUnweighted(G)
				testGraphs(G, G1)

	def testAppend(self):
		n1, n2 = 100, 50
		p1, p2 = 0.01, 0.05
		nodesToDelete = 20

		def testGraphs(G, G1, G2):
			self.assertEqual(G.numberOfNodes(), G1.numberOfNodes() + G2.numberOfNodes())
			self.assertEqual(G.numberOfEdges(), G1.numberOfEdges() + G2.numberOfEdges())
			self.assertEqual(G.isDirected(), G1.isDirected())
			self.assertEqual(G.isDirected(), G2.isDirected())
			self.assertEqual(G.isWeighted(), G1.isWeighted())
			self.assertEqual(G.isWeighted(), G2.isWeighted())

			nodeMap = {}
			v = G1.upperNodeIdBound()
			for u in range(G2.upperNodeIdBound()):
				if G2.hasNode(u):
					nodeMap[u] = v
					v += 1

			G1.forNodes(lambda u: self.assertTrue(G.hasNode(u)))
			G1.forEdges(lambda u, v, w, eid: self.assertTrue(G.hasEdge(u, v)))
			G2.forNodes(lambda u: self.assertTrue(G.hasNode(nodeMap[u])))
			G2.forEdges(lambda u, v, w, eid: self.assertTrue(G.hasEdge(nodeMap[u], nodeMap[v])))

		for seed in range(1, 4):
			nk.setSeed(seed, False)
			random.seed(seed)
			for directed in [True, False]:
				for weighted in [True, False]:
					G1 = nk.generators.ErdosRenyiGenerator(n1, p1, directed).generate()
					G2 = nk.generators.ErdosRenyiGenerator(n2, p2, directed).generate()
					if weighted:
						G1 = self.generateRandomWeights(G1)
						G2 = self.generateRandomWeights(G2)

					G = copy(G1)
					nk.graphtools.append(G, G2)
					testGraphs(G, G1, G2)

					for _ in range(nodesToDelete):
						G1.removeNode(G1.randomNode())
						G2.removeNode(G2.randomNode())
						G3 = copy(G1)
						nk.graphtools.append(G3, G2)
						testGraphs(G3, G1, G2)

	def testMerge(self):
		n1, n2 = 100, 150
		p1, p2 = 0.01, 0.05

		def testGraphs (Gorig, Gmerge, G1):
			for u in range(max(Gorig.upperNodeIdBound(), G1.upperNodeIdBound())):
				self.assertEqual(Gmerge.hasNode(u), Gorig.hasNode(u) or G1.hasNode(u))

			Gorig.forEdges(lambda u, v, w, eid: self.assertTrue(Gmerge.hasEdge(u, v)))
			G1.forEdges(lambda u, v, w, eid: self.assertTrue(Gmerge.hasEdge(u, v)))

			def checkEdges(u, v, w, eid):
				if Gorig.hasNode(u) and Gorig.hasNode(v) and Gorig.hasEdge(u, v):
					self.assertEqual(Gorig.weight(u, v), w)
				else:
					self.assertEqual(G1.weight(u, v), w)
			Gmerge.forEdges(checkEdges)

		for seed in range(1, 4):
			nk.setSeed(seed, False)
			random.seed(seed)
			for directed in [True, False]:
				for weighted in [True, False]:
					Gorig = nk.generators.ErdosRenyiGenerator(n1, p1, directed).generate()
					G1 = nk.generators.ErdosRenyiGenerator(n2, p2, directed).generate()
					if weighted:
						Gorig = self.generateRandomWeights(Gorig)
						G1 = self.generateRandomWeights(G1)
					Gmerge = copy(Gorig)
					nk.graphtools.merge(Gmerge, G1)
					testGraphs(Gorig, Gmerge, G1)

if __name__ == "__main__":
	unittest.main()
