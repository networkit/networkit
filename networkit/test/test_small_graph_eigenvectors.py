import unittest

import networkit as nk
import numpy as np


class TestSmallGraphEigenvectors(unittest.TestCase):
    def testSmallGraphSpectra(self):
        for n in (1, 2):
            for directed in (False, True):
                graph = nk.Graph(n, weighted=True, directed=directed)
                if n == 1:
                    graph.addEdge(0, 0, 3.0)
                else:
                    graph.addEdge(0, 1, 2.0)
                    if directed:
                        graph.addEdge(1, 0, 0.5)
                for kind in ("adjacency", "laplacian"):
                    matrix = getattr(nk.algebraic, kind + "Matrix")(graph).toarray()
                    reference = np.sort(np.linalg.eigvals(matrix))
                    for reverse in (False, True):
                        for cutoff in (-1, 0, n - 1):
                            with self.subTest(n=n, directed=directed, kind=kind, reverse=reverse, cutoff=cutoff):
                                values, vectors = getattr(nk.algebraic, kind + "Eigenvectors")(
                                    graph, cutoff=cutoff, reverse=reverse)
                                count = 1 if cutoff == -1 else cutoff + 1
                                expected = reference[:count] if reverse else reference[-count:]
                                np.testing.assert_allclose(values, expected, atol=1e-12)
                                for value, vector in zip(values, vectors):
                                    np.testing.assert_allclose(matrix @ vector, value * vector, atol=1e-12)

    def testSmallDenseMatrices(self):
        for matrix in (np.array([[2.0]]), np.array([[2.0, 1.0], [1.0, 3.0]])):
            for func in (nk.algebraic.symmetricEigenvectors, nk.algebraic.eigenvectors):
                values, vectors = func(matrix)
                np.testing.assert_allclose(values, [np.linalg.eigvalsh(matrix)[-1]])
                np.testing.assert_allclose(matrix @ vectors[0], values[0] * vectors[0], atol=1e-12)

    def testSingletonCentrality(self):
        centrality = nk.centrality.SciPyEVZ(nk.Graph(1), normalized=True).run()
        np.testing.assert_allclose(centrality.scores(), [1.0])
