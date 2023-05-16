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

	def test_Matrices(self):
		G = nk.readGraph("input/jazz2_directed.gml",nk.Format.GML)
		A = nk.algebraic.adjacencyMatrix(G, "dense")
		self.assertEqual(A[0][0], 1.0)
		self.assertEqual(A[1][1], 1.0)
		self.assertEqual(A[2][2], 0.0)
		self.assertEqual(A[3][3], 0.0)
		self.assertEqual(A[4][4], 0.0)

		B = nk.algebraic.laplacianMatrix(G)
		#B noch asserten

		eigenB = nk.algebraic.laplacianEigenvectors(G)
		self.assertAlmostEqual(eigenB[0][0], -2.2979478724907566e-18-9.534970958477811e-18j, 5)		

	def test_eigenvectorValues(self):
		A_val = np.array([ [-1.3, 2.7, 0.2], [0.8, 4.1, 2.2], [2.1, 4.4, -1.9] ])
		A_mat = scipy.sparse.csr_matrix(A_val)
		eigen = nk.algebraic.eigenvectors(A_mat)
		self.assertAlmostEqual(eigen[0][0], 5.891378730028431+0j, 3)


if __name__ == "__main__":
	unittest.main()
