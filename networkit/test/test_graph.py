#!/usr/bin/env python3
import unittest
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

	def testSubgraphFromNodes(self):
		# Directed
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

		# Undirected
		G = G.toUndirected()

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


if __name__ == "__main__":
	unittest.main()
