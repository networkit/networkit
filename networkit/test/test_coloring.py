#!/usr/bin/env python3
import unittest
import os
import networkit as nk

class TestColoring(unittest.TestCase):

	def testSpectralColoring(self):
		G = nk.readGraph("input/karate.graph", nk.Format.METIS)
		spCol = nk.coloring.SpectralColoring(G)

		spCol.run()

		self.assertLessEqual(len(spCol.getColoring()), G.upperNodeIdBound())

	def checkProperColoring(self, G):
		spCol = nk.coloring.SpectralColoring(G)
		spCol.run()
		coloring = spCol.getColoring()
		self.assertEqual(set(coloring), set(G.iterNodes()))
		for u, v in G.iterEdges():
			self.assertNotEqual(coloring[u], coloring[v])

	def testSpectralColoringSmallGraphs(self):
		for n in [2, 4]:
			G = nk.Graph(n)
			for i in range(n if n > 2 else 1):
				G.addEdge(i, (i + 1) % n)
			self.checkProperColoring(G)

		G = nk.Graph(4, weighted=True)
		for i in range(4):
			G.addEdge(i, (i + 1) % 4, 2.5)
		self.checkProperColoring(G)

	def testSpectralColoringDirected(self):
		G = nk.Graph(20, directed=True)
		for i in range(20):
			G.addEdge(i, (i + 1) % 20)
			G.addEdge(i, (i + 7) % 20)
		self.checkProperColoring(G)

	def testSpectralColoringExhaustedEigenvectors(self):
		spCol = nk.coloring.SpectralColoring(nk.Graph(2))
		spCol.eigenvectors = []
		spCol.colors = {0: set([0, 1])}
		spCol.nextColor = 1
		with self.assertRaises(RuntimeError):
			spCol.split(0)

if __name__ == "__main__":
	unittest.main()
