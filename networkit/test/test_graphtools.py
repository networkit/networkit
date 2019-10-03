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

if __name__ == "__main__":
	unittest.main()
