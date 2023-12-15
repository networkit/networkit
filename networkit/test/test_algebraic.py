#!/usr/bin/env python3
import unittest
import os
import networkit as nk
import numpy as np
import scipy
import random

class TestAlgebraic(unittest.TestCase):

	def testAdjacencyMatrixDenseDirected(self):
		G = nk.Graph(2, directed=True)
		G.addEdge(0,1)
		A = nk.algebraic.adjacencyMatrix(G, "dense")
		self.assertEqual(A[0][0], 0.0)
		self.assertEqual(A[0][1], 1.0)
		self.assertEqual(A[1][0], 0.0)
		self.assertEqual(A[1][1], 0.0)
		self.assertIsInstance(A, np.ndarray)

		G2 = nk.Graph(2, directed=True, weighted=False)
		G2.addEdge(0,1)
		A2 = nk.algebraic.adjacencyMatrix(G2, "dense")
		self.assertEqual(A[0][0], 0)
		self.assertEqual(A[0][1], 1)
		self.assertEqual(A[1][0], 0)
		self.assertEqual(A[1][1], 0)

	def testAdjacencyMatrixDenseUndirected(self):
		G = nk.Graph(2, directed=False)
		G.addEdge(0,1)
		A = nk.algebraic.adjacencyMatrix(G, "dense")
		self.assertEqual(A[0][0], 0.0)
		self.assertEqual(A[0][1], 1.0)
		self.assertEqual(A[1][0], 1.0)
		self.assertEqual(A[1][1], 0.0)
		self.assertIsInstance(A, np.ndarray)

		G2 = nk.Graph(2, directed=False, weighted=False)
		G2.addEdge(0,1)
		A2 = nk.algebraic.adjacencyMatrix(G2, "dense")
		self.assertEqual(A[0][0], 0)
		self.assertEqual(A[0][1], 1)
		self.assertEqual(A[1][0], 1)
		self.assertEqual(A[1][1], 0)
		
	def testAdjacencyMatrixSparse(self):	
		G = nk.Graph(3)
		G.addEdge(0,1)
		G.addEdge(1,2)
		A2 = nk.algebraic.adjacencyMatrix(G, "sparse")
		self.assertIsInstance(A2, scipy.sparse.csr_matrix)

	def testAdjacencyEigenvector(self):
		G = nk.readGraph("input/jazz2_directed.gml",nk.Format.GML)
		eigen1 = nk.algebraic.adjacencyEigenvector(G, 1)
		self.assertAlmostEqual(eigen1[0], 1.0000000000000004+7.45058059692383e-09j, delta=0.5)

	def testLaplacianMatrix(self):	
		G = nk.readGraph("input/jazz2_directed.gml",nk.Format.GML)
		B = nk.algebraic.laplacianMatrix(G)
		self.assertIsInstance(B, scipy.sparse.coo_matrix)

	def testLaplacianEigenvectorsDirected(self):	
		G = nk.readGraph("input/jazz2_directed.gml",nk.Format.GML)
		eigen = nk.algebraic.laplacianEigenvectors(G)
		self.assertAlmostEqual(eigen[0][0], 3.4994785119452635e-34+0j, delta=0.5)
		eigenIndex1 = nk.algebraic.laplacianEigenvector(G, 1)
		self.assertAlmostEqual(eigenIndex1[0], 8.15657495e-01+0.j, delta=0.5)
		
	def testLaplacianEigenvectorsUndirected(self):
		G2 = nk.readGraph("input/jazz2_undirected.gml",nk.Format.GML)
		eigen = nk.algebraic.laplacianEigenvectors(G2)
		self.assertAlmostEqual(eigen[0][0], 8.071939453165987e-18, delta=0.5)
		eigenIndex1 = nk.algebraic.laplacianEigenvector(G2, 1)
		self.assertAlmostEqual(eigenIndex1[0], 3.0, delta=0.5)

	def testEigenvectors(self):
		A = scipy.sparse.csr_matrix(np.array([[-1.3, 2.7, 0.2], [0.8, 4.1, 2.2], [2.1, 4.4, -1.9]]))
		eigen = nk.algebraic.eigenvectors(A)
		self.assertAlmostEqual(eigen[0][0], 5.89137873002843+0j, delta=0.5)

	def testEigenvectorsReverse(self):
		A = scipy.sparse.csr_matrix(np.array([[-1.3, 2.7, 0.2], [0.8, 4.1, 2.2], [2.1, 4.4, -1.9]]))
		eigen = nk.algebraic.eigenvectors(A, reverse=True)
		self.assertAlmostEqual(eigen[0][0], -2.495689365014214+0.517336521966834j, delta=0.5)	

if __name__ == "__main__":
	unittest.main()