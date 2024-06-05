#!/usr/bin/env python3

import networkit as nk
from networkit import vizbridges
import unittest
from collections import defaultdict


class TestVizbridges(unittest.TestCase):
    def getSmallGraph(self, weighted=False, directed=False):
        G = nk.Graph(4, weighted, directed)
        G.addEdge(0, 1, 1.0)
        G.addEdge(0, 2, 2.0)
        G.addEdge(3, 1, 4.0)
        G.addEdge(3, 2, 5.0)
        G.addEdge(1, 2, 3.0)
        if directed:
            G.addEdge(2, 1, 6.0)

        return G

    def testVizNodeScores(self):
        for dim, directed, weighted in zip(
            vizbridges.Dimension, [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted):
                G = self.getSmallGraph(weighted, directed)
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    nodeScores=list(range(G.numberOfNodes())),
                )

    def testVizNodePartition(self):
        for dim, directed, weighted in zip(
            vizbridges.Dimension, [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted):
                G = self.getSmallGraph(weighted, directed)               
                partition = nk.Partition(G.numberOfNodes())
                partition.allToSingletons()
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    nodePartition=partition,
                )

    def testVizNodePalette(self):
        for dim, directed, weighted in zip(
            vizbridges.Dimension, [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted):
                G = self.getSmallGraph(weighted, directed)
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    nodeScores=list(range(G.numberOfNodes())),
                    nodePalette=[(i / 10, i / 10, i / 10) for i in range(4)],
                )

    def testVizShowIds(self):
        for dim, directed, weighted, show in zip(
            vizbridges.Dimension, [True, False], [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted, show=show):
                G = self.getSmallGraph(weighted, directed)
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    nodeScores=list(range(G.numberOfNodes())),
                    showIds=show,
                )

    def testVizCustomSize(self):
        for dim, directed, weighted in zip(
            vizbridges.Dimension, [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted):
                G = self.getSmallGraph(weighted, directed)
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    nodeScores=list(range(G.numberOfNodes())),
                    customSize=50,
                )

    def testVizEdgeScoresList(self):
        for dim, directed, weighted in zip(
            vizbridges.Dimension, [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted):
                G = self.getSmallGraph(weighted, directed)
                G.indexEdges()
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    edgeScores=list(range(G.numberOfEdges())),
                )

    def testVizEdgeScoresDict(self):
        for dim, directed, weighted in zip(
            vizbridges.Dimension, [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted):
                G = self.getSmallGraph(weighted, directed)
                G.indexEdges()
                scores = {}
                for i, [u, v] in enumerate(G.iterEdges()):
                    scores[u, v] = i
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    edgeScores=scores,
                )

    def testVizEdgeScoresUndirectedReversed(self):
        for dim, weighted in zip(
            vizbridges.Dimension, [True, False]
        ):
            with self.subTest(dim=dim, weighted=weighted):
                G = self.getSmallGraph(weighted, False)
                G.indexEdges()
                scores = {}
                for i, [u, v] in enumerate(G.iterEdges()):
                    scores[v, u] = i
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    edgeScores=scores,
                )

    def testVizEdgeScoresDefaultdict(self):
        for dim, directed, weighted in zip(
            vizbridges.Dimension, [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted):
                G = self.getSmallGraph(weighted, directed)
                G.indexEdges()
                scores = defaultdict(lambda: 0)
                for u, v in G.iterEdges():
                    scores[u, v] = 1
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    edgeScores=scores,
                )

    def testVizEdgeScoresDefaultdictOneEntry(self):
        for dim, directed, weighted in zip(
            vizbridges.Dimension, [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted):
                G = self.getSmallGraph(weighted, directed)
                G.indexEdges()
                scores = defaultdict(lambda: 0)
                # create incomplete defaultdict
                for u, v in G.iterEdges():
                    scores[u, v] = 1
                    break
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    edgeScores=scores,
                )

    def testVizEdgePalette(self):
        for dim, directed, weighted in zip(
            vizbridges.Dimension, [True, False], [True, False]
        ):
            with self.subTest(dim=dim, directed=directed, weighted=weighted):
                G = self.getSmallGraph(weighted, directed)
                G.indexEdges()
                vizbridges.widgetFromGraph(
                    G,
                    dimension=dim,
                    edgeScores=list(range(G.numberOfEdges())),
                    edgePalette=[(i / 10, i / 10, i / 10) for i in range(4)],
                )

    def testNodeScoreAndPartitionExclusive(self):
        G = self.getSmallGraph(False, False)
        partition = nk.community.ClusteringGenerator(G).makeRandomClustering(G, 3)
        with self.assertRaises(Exception):
            vizbridges.widgetFromGraph(
                G, nodeScores=list(range(G.numberOfNodes())), nodePalette=partition
            )

    def testCompleteNodeScores(self):
        G = self.getSmallGraph(False, False)
        with self.assertRaises(Exception):
            vizbridges.widgetFromGraph(G, nodeScores=list(range(G.numberOfNodes() - 1)))

    def testCompleteEdgeScores(self):
        G = self.getSmallGraph(False, False)
        with self.assertRaises(Exception):
            vizbridges.widgetFromGraph(G, edgeScores=list(range(G.numberOfEdges() - 1)))

    def testNoneInEdgeScores(self):
        G = self.getSmallGraph(False, False)
        scores = [None] * G.numberOfEdges()
        with self.assertRaises(Exception):
            vizbridges.widgetFromGraph(G, edgeScores=scores)


if __name__ == "__main__":
    unittest.main()
