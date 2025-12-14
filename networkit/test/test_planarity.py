#!/usr/bin/env python3
import unittest

import networkit as nk

class TestLeftRightPlanarityCheck(unittest.TestCase):

    def test_no_edges_indexed_graph_throws(self):
        # C++: Graph graph(0);  (no edge IDs)
        G = nk.Graph(0)  # default: undirected, unweighted; edge IDs typically OFF
        with self.assertRaises(RuntimeError) as ctx:
            nk.planarity.LeftRightPlanarityCheck(G)
        self.assertEqual(str(ctx.exception), "The graph has no edge IDs.")

    def test_directed_graph_throws(self):
        # C++: Graph graph(0, false, true, false); -> directed graph
        # Python: directed=True, and edge ids OFF here doesn't matter because ctor should reject directed anyway.
        G = nk.Graph(0, directed=True)
        with self.assertRaises(RuntimeError) as ctx:
            nk.planarity.LeftRightPlanarityCheck(G)
        self.assertEqual(str(ctx.exception), "The graph is not an undirected graph.")

    def test_is_planar_throws_if_run_not_called(self):
        # C++: Graph graph(0, false, false, true); -> edge IDs enabled
        # Python: edge id indexing is either via indexEdges() or via constructor flag depending on API version.
        # Use indexEdges() explicitly to be safe and to mirror the C++ test setup.
        G = nk.Graph(0)
        G.indexEdges()

        algo = nk.planarity.LeftRightPlanarityCheck(G)
        with self.assertRaises(RuntimeError) as ctx:
            algo.isPlanar()
        self.assertEqual(str(ctx.exception), "Error, run must be called first")

    def test_planar_empty_graph(self):
        G = nk.Graph(0)
        G.indexEdges()
        algo = nk.planarity.LeftRightPlanarityCheck(G).run()
        self.assertTrue(algo.isPlanar())

    def test_planar_single_node(self):
        G = nk.Graph(1)
        G.indexEdges()
        algo = nk.planarity.LeftRightPlanarityCheck(G).run()
        self.assertTrue(algo.isPlanar())

    def test_non_planar_k33(self):
        G = nk.Graph(6)
        edges = [
            (0, 3), (0, 4), (0, 5),
            (1, 3), (1, 4), (1, 5),
            (2, 3), (2, 4), (2, 5),
        ]
        for u, v in edges:
            G.addEdge(u, v)
        G.indexEdges()

        algo = nk.planarity.LeftRightPlanarityCheck(G).run()
        self.assertFalse(algo.isPlanar())

    def test_non_planar_k333(self):
        G = nk.Graph(9)

        # Partitions: A={0,1,2}, B={3,4,5}, C={6,7,8}
        A = [0, 1, 2]
        B = [3, 4, 5]
        C = [6, 7, 8]

        # Complete tripartite: add all edges between different partitions
        for u in A:
            for v in B:
                G.addEdge(u, v)
            for v in C:
                G.addEdge(u, v)
        for u in B:
            for v in C:
                G.addEdge(u, v)

        G.indexEdges()
        algo = nk.planarity.LeftRightPlanarityCheck(G).run()
        self.assertFalse(algo.isPlanar())

    def test_one_planar_one_non_planar_subgraph(self):
        G = nk.Graph(10)

        # K3,3 on {0,1,2} x {3,4,5}
        for u in [0, 1, 2]:
            for v in [3, 4, 5]:
                G.addEdge(u, v)

        # Cycle on {6,7,8,9}
        cycle_edges = [(6, 7), (7, 8), (8, 9), (9, 6)]
        for u, v in cycle_edges:
            G.addEdge(u, v)

        G.indexEdges()
        algo = nk.planarity.LeftRightPlanarityCheck(G).run()
        self.assertFalse(algo.isPlanar())

    def test_planar_4elt_graph(self):
        # This depends on whether METISGraphReader is exposed in Python as nk.readGraph / nk.graphio / nk.METISGraphReader.
        # Common NetworKit Python API: nk.readGraph(path, nk.Format.METIS)
        G = nk.readGraph("input/4elt.graph", nk.Format.METIS)
        G.indexEdges()
        algo = nk.planarity.LeftRightPlanarityCheck(G).run()
        self.assertTrue(algo.isPlanar())


if __name__ == "__main__":
    unittest.main()
