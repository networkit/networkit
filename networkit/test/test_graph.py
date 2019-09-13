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

	def testSubgraphFromNodesDirected(self):
		G = self.getSmallGraph(True, True)

		res = G.subgraphFromNodes([0])
		self.assertTrue(res.isWeighted())
		self.assertTrue(res.isDirected())
		self.assertEqual(res.numberOfNodes(), 1)
		self.assertEqual(res.numberOfEdges(), 0)

		res = G.subgraphFromNodes([0], True)
		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 2)

		res = G.subgraphFromNodes([0, 1])
		self.assertEqual(res.numberOfNodes(), 2)
		self.assertEqual(res.numberOfEdges(), 1)

		res = G.subgraphFromNodes([0, 1], True)
		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 3)

		res = G.subgraphFromNodes([0, 1], True, True)
		self.assertEqual(res.numberOfNodes(), 4)
		self.assertEqual(res.numberOfEdges(), 4)

	def testSubgraphFromNodesUndirected(self):
		G = self.getSmallGraph(True, False)

		res = G.subgraphFromNodes([0])
		self.assertTrue(res.isWeighted())
		self.assertFalse(res.isDirected())
		self.assertEqual(res.numberOfNodes(), 1)
		self.assertEqual(res.numberOfEdges(), 0)

		res = G.subgraphFromNodes([0], True)
		self.assertEqual(res.numberOfNodes(), 3)
		self.assertEqual(res.numberOfEdges(), 2)
		self.assertEqual(G.weight(0, 1), 1.0)
		self.assertEqual(G.weight(0, 2), 2.0)

		res = G.subgraphFromNodes([0, 1])
		self.assertEqual(res.numberOfNodes(), 2)
		self.assertEqual(res.numberOfEdges(), 1)

		res = G.subgraphFromNodes([0, 1], True)
		self.assertEqual(res.numberOfNodes(), 4)
		self.assertEqual(res.numberOfEdges(), 4)

		res = G.subgraphFromNodes(set([0, 1]), True)
		self.assertEqual(res.numberOfNodes(), 4)
		self.assertEqual(res.numberOfEdges(), 4)

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

	def testGraphTranspose(self):
		nk.setSeed(1, True)
		G = nk.generators.ErdosRenyiGenerator(100, 0.2, True).generate()

		for i in range(20):
			u = G.randomNode()
			if not G.hasEdge(u, u):
				G.addEdge(u, u)
		self.assertGreater(G.numberOfSelfLoops(), 0)

		# Delete a few nodes
		for i in range(10):
			G.removeNode(G.randomNode())
		self.assertGreater(G.numberOfSelfLoops(), 0)

		# Assign random weights
		GWeighted = nk.Graph(G, True, True)
		for u, v in GWeighted.edges():
			GWeighted.setWeight(u, v, random.random())

		GWeighted.indexEdges()

		GTrans = GWeighted.transpose()

		for u, v in GWeighted.edges():
			self.assertEqual(GWeighted.edgeId(u, v), GTrans.edgeId(v, u))
			self.assertEqual(GWeighted.weight(u, v), GTrans.weight(v, u))

		for v, u in GTrans.edges():
			self.assertEqual(GWeighted.edgeId(u, v), GTrans.edgeId(v, u))
			self.assertEqual(GWeighted.weight(u, v), GTrans.weight(v, u))

		for u in range(GWeighted.upperNodeIdBound()):
			self.assertEqual(GWeighted.hasNode(u), GTrans.hasNode(u))

		self.assertEqual(GWeighted.upperEdgeIdBound(),  GTrans.upperEdgeIdBound())
		self.assertEqual(GWeighted.upperNodeIdBound(),  GTrans.upperNodeIdBound())
		self.assertEqual(GWeighted.numberOfSelfLoops(), GTrans.numberOfSelfLoops())
		self.assertEqual(GWeighted.numberOfNodes(),     GTrans.numberOfNodes())
		self.assertEqual(GWeighted.numberOfEdges(),     GTrans.numberOfEdges())


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
					print(nodeList)
			except StopIteration:
				pass
			return nodeList

		for node in range(G.upperNodeIdBound()):
			forEdgesNbrs = []
			G.forInEdgesOf(node, nbrFunc)
			nodeInNbrs = nodeInIter(node)
			self.assertEqual(sorted(forEdgesNbrs), sorted(nodeInNbrs))

if __name__ == "__main__":
	unittest.main()
