#!/usr/bin/env python3
import unittest
import os
import networkit as nk
import numpy as np
import scipy


class Test_Centrality(unittest.TestCase):

	def test_DegreeCentrality(self):
		g = nk.Graph(8, False, False)

		g.addEdge(0, 2)
		g.addEdge(0, 5)
		g.addEdge(1, 2)
		g.addEdge(2, 3)
		g.addEdge(2, 2)
		g.addEdge(2, 4)
		g.addEdge(3, 5)
		g.addEdge(4, 5)
		g.addEdge(5, 5)
		g.addEdge(5, 6)
		g.addEdge(5, 7)
		g.addEdge(7, 7)

		expected_result = [2.0, 1.0, 4.0, 2.0, 2.0, 5.0, 1.0, 1.0]

		dc = nk.centrality.DegreeCentrality(g).run().scores()

		self.assertListEqual(expected_result, dc)

	def test_SquareClusteringCoefficient(self):
		g = nk.Graph()
		g.addNodes(7)
		edges = [(0, 1), (1, 2), (2, 3), (0, 3), (3, 4), (4, 5), (5, 6), (6, 3)]
		[g.addEdge(*edge) for edge in edges]

		expected_result = [1 / 3, 1.0, 1 / 3, 0.2, 1 / 3, 1.0, 1 / 3]
		scores = nk.centrality.LocalSquareClusteringCoefficient(g).run().scores()
		self.assertListEqual(expected_result, scores)

	def test_ForestCentralty(self):
		nk.engineering.setSeed(42, False)
		eps = 0.05
		g = nk.generators.HyperbolicGenerator(200).generate()
		root = nk.graphtools.augmentGraph(g)

		fc = nk.centrality.ForestCentrality(g, root, eps)
		fc.run()
		apxDiag = fc.getDiagonal()

		A = nk.algebraic.adjacencyMatrix(g, "dense")
		Fmat = scipy.sparse.csgraph.laplacian(A)
		diag = np.linalg.pinv(Fmat).diagonal()

		for apx, exact in zip(apxDiag, diag):
			self.assertLessEqual(abs(apx - exact), eps)

if __name__ == "__main__":
	unittest.main()
