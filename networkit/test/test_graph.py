#!/usr/bin/env python3
import unittest
import random
import networkit as nk

class TestGraph(unittest.TestCase):
	def getSmallGraph(self, weighted=False, directed=False):
		G = nk.Graph(4, weighted, directed)
		G.addEdge(0, 1, 1.0)
		G.addEdge(0, 2, 2.0)
		G.addEdge(3, 1, 4.0)
		G.addEdge(3, 2, 5.0)
		G.addEdge(1, 2, 3.0)

		return G

	def testRemoveMultiEdges(self):

		def addMultiEdges(G, nMultiEdges):
			for i in range(nMultiEdges):
				e = G.randomEdge()
				G.addEdge(e[0], e[1])
			return G

		def testGraph(G):
			addMultiEdges(G, nMultiEdges)
			self.assertEqual(G.numberOfEdges(), m + nMultiEdges)
			G.removeMultiEdges()
			self.assertEqual(G.numberOfEdges(), m)

		nMultiEdges = 5

		# Directed
		G = self.getSmallGraph(True, True)
		m = G.numberOfEdges()
		testGraph(G)

		# Undirected
		G = self.getSmallGraph(True, False)
		m = G.numberOfEdges()
		testGraph(G)

	def testNeighbors(self):
		# Directed
		G = nk.Graph(4, False, True)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)
		G.addEdge(1, 2)

		self.assertEqual(sorted(G.neighbors(0)), [1, 2])
		self.assertEqual(sorted(G.neighbors(1)), [2])
		self.assertEqual(sorted(G.neighbors(2)), [])
		self.assertEqual(sorted(G.neighbors(3)), [1, 2])

		# Undirected
		G = nk.Graph(4, False, False)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)
		G.addEdge(1, 2)

		self.assertEqual(sorted(G.neighbors(0)), [1, 2])
		self.assertEqual(sorted(G.neighbors(1)), [0, 2, 3])
		self.assertEqual(sorted(G.neighbors(2)), [0, 1, 3])
		self.assertEqual(sorted(G.neighbors(3)), [1, 2])

	def testRandomEdgesReproducibility(self):
		numSamples = 10
		numSeeds = 3
		numRepeats = 10

		for directed in [False, True]:
			G = self.getSmallGraph(False, directed)

			results = [[] for i in range(numSeeds)]
			for repeats in range(numRepeats):
				for seed in range(numSeeds):
					nk.setSeed(seed, False)
					results[seed].append(G.randomEdges(numSamples))

			# assert results are different for different seeds
			for seed in range(1, numSeeds):
				self.assertNotEqual(results[0][0], results[seed][0])

			# assert results are consistent for same seeds
			for repeats in results:
				for x in repeats[1:]:
					self.assertEqual(repeats[0], x)

	def testIterator(self):
		# Undirected
		G = nk.Graph(4, False, False)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(1, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)

		def nbrFunc(u, v, weight, edgeId):
			forEdgesNbrs.append(v)

		# Iterate through neighbours of node using iterNeighbors
		def nodeIter(node):
			nodeList = []
			nbrs = G.iterNeighbors(node)
			try:
				while nbrs is not None:
					nodeList.append(next(nbrs))
			except StopIteration:
				pass
			return nodeList

		for node in range(G.upperNodeIdBound()):
			forEdgesNbrs = []
			G.forEdgesOf(node, nbrFunc)
			nodeNbrs = nodeIter(node)
			self.assertEqual(sorted(forEdgesNbrs), sorted(nodeNbrs))

		# Directed
		G = nk.Graph(4, False, True)

		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(3, 1)
		G.addEdge(3, 2)
		G.addEdge(1, 2)

		# Iterate through neighbours of node using iterNeighbors
		def nodeInIter(node):
			nodeList = []
			nbrs = G.iterInNeighbors(node)
			try:
				while nbrs is not None:
					nodeList.append(next(nbrs))
			except StopIteration:
				pass
			return nodeList

		for node in range(G.upperNodeIdBound()):
			forEdgesNbrs = []
			G.forInEdgesOf(node, nbrFunc)
			nodeInNbrs = nodeInIter(node)
			self.assertEqual(sorted(forEdgesNbrs), sorted(nodeInNbrs))

	def testAddNodes(self):
		G = nk.Graph(0)
		G.addNodes(10)
		self.assertEqual(G.numberOfNodes(), 10)
		for u in range(G.numberOfNodes()):
			self.assertTrue(G.hasNode(u))

	def testAddEdgeWithMissingNodes(self):
		G = nk.Graph(0)

		with self.assertRaises(RuntimeError):
			G.addEdge(0, 1)

		self.assertEqual(G.numberOfNodes(), 0)
		self.assertEqual(G.numberOfEdges(), 0)

		G.addEdge(0, 2, addMissing = True)

		self.assertEqual(G.numberOfNodes(), 3)
		self.assertEqual(G.numberOfEdges(), 1)

		G.removeNode(1)

		self.assertEqual(G.numberOfNodes(), 2)
		self.assertEqual(G.numberOfEdges(), 1)

		G.addEdge(0, 1, addMissing = True)

		self.assertEqual(G.numberOfNodes(), 3)
		self.assertEqual(G.numberOfEdges(), 2)

		G.removeNode(2)

		G.addEdge(1, 2, addMissing = True)

		self.assertEqual(G.numberOfNodes(), 3)
		self.assertEqual(G.numberOfEdges(), 2)

	def testMaxDegree(self):
		n = 100
		p = 0.2
		edgeUpdates = 10

		def computeMaxDeg(G, inDegree = False):
			nodes = []
			G.forNodes(lambda u: nodes.append(u))
			maxDeg = 0
			for u in nodes:
				maxDeg = max(maxDeg, G.degreeIn(u) if inDegree else G.degreeOut(u))
			return maxDeg

		def doTest(G):
			self.assertEqual(G.maxDegree(), computeMaxDeg(G))
			self.assertEqual(G.maxDegreeIn(), computeMaxDeg(G, True))

		for seed in range(1, 4):
			nk.setSeed(seed, False)
			for directed in [True, False]:
				for weighted in [True, False]:
					G = nk.generators.ErdosRenyiGenerator(n, p, directed).generate()
					if weighted:
						G = nk.graphtools.toWeighted(G)

					doTest(G)
					for _ in range(edgeUpdates):
						e = G.randomEdge()
						G.removeEdge(e[0], e[1])
						doTest(G)

					for _ in range(edgeUpdates):
						e = G.randomNode(), G.randomNode()
						while G.hasEdge(e[0], e[1]):
							e = G.randomNode(), G.randomNode()
						G.addEdge(e[0], e[1])
						doTest(G)

	def testWeightedDegree(self):
		n = 100
		p = 0.2

		for seed in range(1, 4):
			nk.setSeed(seed, False)
			random.seed(seed)
			for directed in [True, False]:
				for weighted in [True, False]:
					G = nk.generators.ErdosRenyiGenerator(n, p, directed).generate()
					if weighted:
						G = nk.graphtools.toWeighted(G)
						G.forEdges(lambda u, v, w, eid: G.setWeight(u, v, random.random()))

					def testWeightedDegreeOfNode(u):
						wDeg, wDegTwice = 0, 0
						for v in G.iterNeighbors(u):
							w = G.weight(u, v)
							wDeg += w
							wDegTwice += w if u != v else 2 * w

						self.assertEqual(G.weightedDegre(u), wDeg)
						self.assertEqual(G.weightedDegre(u, True), wDegTwice)

						wInDeg, wInDegTwice = 0, 0
						for v in G.iterInNeighbors(u):
							w = G.weight(v, u)
							wInDeg += w
							wInDegTwice += w if u != v else 2 * w

						self.assertEqual(G.weightedDegreeIn(u), wInDeg)
						self.assertEqual(G.weightedDegreeIn(u, True), wInDegTwice)

if __name__ == "__main__":
	unittest.main()
