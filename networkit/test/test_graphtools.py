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
		for e in G.iterEdges():
			G.setWeight(e[0], e[1], random.random())

		return G

	def testMaxDegree(self):
		n = 100
		p = 0.2
		edgeUpdates = 10

		def computeMaxDeg(G, weighted = False, inDegree = False):
			nodes = []
			G.forNodes(lambda u: nodes.append(u))
			maxDeg = 0
			def getDegree(u):
				if weighted:
					return G.weightedDegreeIn(u) if inDegree else G.weightedDegree(u)
				return G.degreeIn(u) if inDegree else G.degreeOut(u)

			for u in nodes:
				maxDeg = max(maxDeg, getDegree(u))

			return maxDeg

		def doTest(G):
			self.assertEqual(nk.graphtools.maxDegree(G), computeMaxDeg(G, False))
			self.assertEqual(nk.graphtools.maxInDegree(G), computeMaxDeg(G, False, True))
			self.assertEqual(nk.graphtools.maxWeightedDegree(G), computeMaxDeg(G, True))
			self.assertEqual(nk.graphtools.maxWeightedInDegree(G), computeMaxDeg(G, True, True))

		for seed in range(1, 4):
			nk.setSeed(seed, False)
			for directed in [True, False]:
				for weighted in [True, False]:
					G = nk.generators.ErdosRenyiGenerator(n, p, directed).generate()
					if weighted:
						G = nk.graphtools.toWeighted(G)

					doTest(G)
					for _ in range(edgeUpdates):
						e = nk.graphtools.randomEdge(G)
						G.removeEdge(e[0], e[1])
						doTest(G)

					for _ in range(edgeUpdates):
						e = nk.graphtools.randomNode(G), nk.graphtools.randomNode(G)
						while G.hasEdge(e[0], e[1]):
							e = nk.graphtools.randomNode(G), nk.graphtools.randomNode(G)
						G.addEdge(e[0], e[1])
						doTest(G)

	def testRandomNode(self):
		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)
				for i in range(10):
					self.assertTrue(G.hasNode(nk.graphtools.randomNode(G)))
				n = G.numberOfNodes()
				for i in range(n):
					G.removeNode(i)
				self.assertEqual(nk.graphtools.randomNode(G), nk.none)

	def testRandomNeighbor(self):
		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)
				for i in range(10):
					u = nk.graphtools.randomNode(G)
					v = nk.graphtools.randomNeighbor(G, u)
					self.assertNotEqual(G.degree(u) == 0, G.hasNode(v))

	def testRandomEdge(self):
		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)
				for i in range(10):
					u, v = nk.graphtools.randomEdge(G)
					self.assertTrue(G.hasEdge(u, v))

	def testRandomEdges(self):
		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)
				for i in range(10):
					randomEdges = nk.graphtools.randomEdges(G, 5)
					for u, v in randomEdges:
						self.assertTrue(G.hasEdge(u, v))

	def testSize(self):
		def doTest(G):
			(n, m) = nk.graphtools.size(G)
			self.assertEqual(G.numberOfNodes(), n)
			self.assertEqual(G.numberOfEdges(), m)

		n, p = 100, 0.1
		for seed in range(1, 4):
			for directed in [True, False]:
				for weighted in [True, False]:
					G = nk.generators.ErdosRenyiGenerator(n, p, directed).generate()
					if weighted:
						G = nk.graphtools.toWeighted(G)
					doTest(G)

					for _ in range(10):
						G.removeNode(nk.graphtools.randomNode(G))
						doTest(G)

					for _ in range(10):
						if G.numberOfEdges() == 0:
							break
						u, v = nk.graphtools.randomEdge(G)
						G.removeEdge(u, v)
						doTest(G)

	def testVolume(self):
		for directed in [True, False]:
			for weighted in [True, False]:
				G = self.getSmallGraph(weighted, directed)

				# Test total volume
				v1 = nk.graphtools.volume(G)
				mod = 1.0 if G.isDirected() else 2.0;
				baseG = G.totalEdgeWeight() if G.isWeighted() else G.numberOfEdges();
				self.assertEqual(v1, mod * baseG)

				# Test partial volume
				nodes = [nk.graphtools.randomNode(G)]
				v2 = nk.graphtools.inVolume(G, nodes)
				v3 = nk.graphtools.volume(G, nodes)
				self.assertLessEqual(v2, v1)
				self.assertLessEqual(v3, v1)

	def testDensity(self):
		def doTest(G):
			d = nk.graphtools.density(G)
			self.assertGreaterEqual(d, 0)
			self.assertEqual(
					d > 0,
					G.numberOfNodes() >= 1 and G.numberOfEdges() - G.numberOfSelfLoops() > 0)

		n, p = 100, 0.1
		for seed in range(1, 4):
			for directed in [True, False]:
				for weighted in [True, False]:
					G = nk.generators.ErdosRenyiGenerator(n, p, directed).generate()
					if weighted:
						G = nk.graphtools.toWeighted(G)
					doTest(G)

					for _ in range(n):
						G.removeNode(nk.graphtools.randomNode(G))
						doTest(G)

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
					G.removeNode(nk.graphtools.randomNode(G))
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

		res = nk.graphtools.subgraphAndNeighborsFromNodes(G, nodes, True)

		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 2) # 0-1, 0-2, NOT 1-2

		self.assertEqual(G.weight(0, 1), 1.0)
		self.assertEqual(G.weight(0, 2), 2.0)

		nodes = set([0, 1])
		res = nk.graphtools.subgraphFromNodes(G, nodes)
		self.assertEqual(res.numberOfNodes(), 2)
		self.assertEqual(res.numberOfEdges(), 1) # 0 - 1

		res = nk.graphtools.subgraphAndNeighborsFromNodes(G, nodes, True)
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
		res = nk.graphtools.subgraphAndNeighborsFromNodes(G, nodes, True)
		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 2) # 0->1, 0->2, NOT 1->2

		nodes = set([0, 1])
		res = nk.graphtools.subgraphFromNodes(G, nodes)
		self.assertEqual(res.numberOfNodes(), 2)
		self.assertEqual(res.numberOfEdges(), 1) # 0 -> 1

		nodes = set([0, 1])
		res = nk.graphtools.subgraphAndNeighborsFromNodes(G, nodes, True)
		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 3) # 0->1, 0->2, 1->2

		nodes = set([0, 1])
		res = nk.graphtools.subgraphAndNeighborsFromNodes(G, nodes, True, True)
		self.assertEqual(res.numberOfNodes(), 4)
		self.assertEqual(res.numberOfEdges(), 4) # 0->1, 0->2, 1->2, 3->1

	def testGraphTranspose(self):
		for seed in range(1, 4):
			nk.setSeed(seed, True)
			random.seed(seed)
			G = nk.generators.ErdosRenyiGenerator(100, 0.2, True).generate()

			for _ in range(20):
				u = nk.graphtools.randomNode(G)
				if not G.hasEdge(u, u):
					G.addEdge(u, u)

			# Delete a few nodes
			for _ in range(10):
				G.removeNode(nk.graphtools.randomNode(G))
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
			edgesLost = 0;
			for e in G.iterEdges():
				if G.hasEdge(e[1], e[0]) and e[1] != e[0]:
					edgesLost += 1
			self.assertEqual(G.numberOfEdges() - edgesLost, G1.numberOfEdges())
			self.assertEqual(G.upperEdgeIdBound(), G1.upperEdgeIdBound())
			self.assertEqual(G.isWeighted(), G1.isWeighted())
			self.assertNotEqual(G.isDirected(), G1.isDirected())
			self.assertEqual(G.hasEdgeIds(), G1.hasEdgeIds())

			def testEdges(u, v, w, eid):
				self.assertTrue(G1.hasEdge(u, v))
				self.assertEqual(G1.weight(u, v), w if not G.hasEdge(v, u) else w + G.weight(v, u))
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
						G1.removeNode(nk.graphtools.randomNode(G1))
						G2.removeNode(nk.graphtools.randomNode(G2))
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

	def testRemoveEdgesFromIsolatedSet(self):
		n = 6

		def generateTwoComponents(directed, weighted):
			G = nk.Graph(n, directed, weighted)
			G.addEdge(0, 1)
			G.addEdge(1, 2)
			G.addEdge(2, 0)

			G.addEdge(3, 4)
			G.addEdge(4, 5)
			G.addEdge(5, 3)

			return G

		for directed in [True, False]:
			for weighted in [True, False]:
				G = generateTwoComponents(directed, weighted)
				nk.graphtools.removeEdgesFromIsolatedSet(G, [0, 1, 2])
				self.assertEqual(G.numberOfEdges(), 3)
				self.assertTrue(G.hasEdge(3, 4))
				self.assertTrue(G.hasEdge(4, 5))
				self.assertTrue(G.hasEdge(5, 3))
				nk.graphtools.removeEdgesFromIsolatedSet(G, [3, 4, 5])
				self.assertEqual(G.numberOfEdges(), 0)

	def testSortEdgesByWeight(self):
		def checkSortedEdges(g, decreasing):
			for u in g.iterNodes():
				prevNode = g.numberOfNodes() if decreasing else -1
				prevWeight = 2 if decreasing else 0

				for v in g.iterNeighbors(u):
					w = g.weight(u, v)
					if decreasing:
						if w == prevWeight:
							self.assertLess(prevNode, v)
						else:
							self.assertLess(w, prevWeight)
					else:
						if w == prevWeight:
							self.assertLess(prevNode, v)
						else:
							self.assertGreater(w, prevWeight)

					prevNode, prevWeight = v, w

		def doTest(g):
			nk.graphtools.sortEdgesByWeight(g, False)
			checkSortedEdges(g, False)
			nk.graphtools.sortEdgesByWeight(g, True)
			checkSortedEdges(g, True)

		g = nk.readGraph('input/PGPgiantcompo.graph', nk.Format.METIS)
		g.removeSelfLoops()
		g.removeMultiEdges()

		# Test unweighted
		doTest(g)

		random.seed(1)
		g = self.generateRandomWeights(g)
		e = nk.graphtools.randomEdge(g)

		# Test weighted
		doTest(g)

	def testTopologicalSort(self):
		def generateGraph(directed):
			G = nk.Graph(5, False, directed)
			G.addEdge(0, 1)
			G.addEdge(0, 2)
			G.addEdge(2, 1)
			G.addEdge(1, 3)
			G.addEdge(4, 2)
			return G
		
		for directed in [True, False]:
			G = generateGraph(directed)
			if(directed == False):
				with self.assertRaises(Exception):
					nk.graphtools.topologicalSort(G)
			else:
				res = nk.graphtools.topologicalSort(G)
				indexNode0 = res.index(0)
				indexNode1 = res.index(1)
				indexNode2 = res.index(2)
				indexNode3 = res.index(3)
				indexNode4 = res.index(4)
				self.assertLess(indexNode0, indexNode2)
				self.assertLess(indexNode4, indexNode2)
				self.assertLess(indexNode2, indexNode1)
				self.assertLess(indexNode1, indexNode3)

	def testAugmentedGraph(self):
		nk.engineering.setSeed(42, False)
		n = 20
		G = nk.generators.ErdosRenyiGenerator(n, 0.3).generate()
		augG, root = nk.graphtools.createAugmentedGraph(G);

		self.assertTrue(augG.hasNode(root));
		self.assertEqual(n + 1, augG.numberOfNodes());
		self.assertEqual(G.numberOfEdges() + n, augG.numberOfEdges());
		self.assertEqual(n, augG.degree(root));

if __name__ == "__main__":
	unittest.main()
