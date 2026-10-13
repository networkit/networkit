import unittest

import networkit as nk
from networkit import vizbridges


class TestLinearScoreColors(unittest.TestCase):
    def testPublicWidgetScoreColors(self):
        graph = nk.Graph(5)
        for i in range(5):
            graph.addEdge(i, (i + 1) % 5)
        graph.indexEdges()
        palette = [(i / 255, 0.0, 0.0) for i in range(256)]
        scores = [-2.0, -1.0, 0.0, 1.0, 2.0]
        widget = vizbridges.widgetFromGraph(
            graph, nodeScores=scores, edgeScores=scores, nodePalette=palette, edgePalette=palette)
        expected = [(0, 0, 0), (63, 0, 0), (127, 0, 0), (191, 0, 0), (255, 0, 0)]
        self.assertEqual([tuple(node.data["color"]) for node in widget.graph.nodes], expected)
        for edge in widget.graph.edges:
            edge_id = graph.edgeId(edge.data["source"], edge.data["target"])
            self.assertEqual(tuple(edge.data["color"]), expected[edge_id])
