#!/usr/bin/env python3
import unittest
import unittest.mock
import io
import networkit as nk
import networkx as nx
import warnings

from typing import Union, List


class TestNXAdapter(unittest.TestCase):
    def toyNxgraph(
        self, directed: bool, data: List[str]
    ) -> Union[nx.Graph, nx.DiGraph]:
        nxGraph = nx.Graph()

        # 0 - 'a'
        # |    |
        # 3    2
        nxGraph.add_edge(0, "a")
        nxGraph.add_edge(0, 3)
        nxGraph.add_edge("a", 2)

        if "weight" in data:
            nxGraph[0]["a"]["weight"] = 2
            nxGraph[0][3]["weight"] = 4
            nxGraph["a"][2]["weight"] = -6

        # int attribute
        if "intEdgeAttr" in data:
            nxGraph[0]["a"]["intEdgeAttr"] = 1
            nxGraph[0][3]["intEdgeAttr"] = 2
            nxGraph["a"][2]["intEdgeAttr"] = -3

        if "intNodeAttr" in data:
            nxGraph.nodes[0]["intNodeAttr"] = 4
            nxGraph.nodes["a"]["intNodeAttr"] = 5
            nxGraph.nodes[2]["intNodeAttr"] = 6
            nxGraph.nodes[3]["intNodeAttr"] = 7

        # float attribute
        if "floatEdgeAttr" in data:
            nxGraph[0]["a"]["floatEdgeAttr"] = 1.2
            nxGraph[0][3]["floatEdgeAttr"] = 2.3
            nxGraph["a"][2]["floatEdgeAttr"] = -3.1

        if "floatNodeAttr" in data:
            nxGraph.nodes[0]["floatNodeAttr"] = 4.7
            nxGraph.nodes["a"]["floatNodeAttr"] = 5.6
            nxGraph.nodes[2]["floatNodeAttr"] = 6.5
            nxGraph.nodes[3]["floatNodeAttr"] = 7.4

        # str attribute
        if "strEdgeAttr" in data:
            nxGraph[0]["a"]["strEdgeAttr"] = "0a"
            nxGraph[0][3]["strEdgeAttr"] = "03"
            nxGraph["a"][2]["strEdgeAttr"] = "a2"

        if "strNodeAttr" in data:
            nxGraph.nodes[0]["strNodeAttr"] = "n0"
            nxGraph.nodes["a"]["strNodeAttr"] = "na"
            nxGraph.nodes[2]["strNodeAttr"] = "n2"
            nxGraph.nodes[3]["strNodeAttr"] = "n3"

        # complex attribute: tuple
        if "complexEdgeAttr" in data:
            nxGraph[0]["a"]["complexEdgeAttr"] = (0, "a", "0a")
            nxGraph[0][3]["complexEdgeAttr"] = (0, 3, "03")
            nxGraph["a"][2]["complexEdgeAttr"] = ("a", 2, "a2")

        if "complexNodeAttr" in data:
            nxGraph.nodes[0]["complexNodeAttr"] = (0, "0")
            nxGraph.nodes["a"]["complexNodeAttr"] = ("a", "a")
            nxGraph.nodes[2]["complexNodeAttr"] = (2, "2")
            nxGraph.nodes[3]["complexNodeAttr"] = (3, "3")

        if directed:
            nxGraph = nxGraph.to_directed()
            nxGraph.add_edge(3, 2)

            if "weight" in data:
                nxGraph[3][2]["weight"] = -5
            if "intEdgeAttr" in data:
                nxGraph[3][2]["intEdgeAttr"] = -1
            if "floatEdgeAttr" in data:
                nxGraph[3][2]["floatEdgeAttr"] = -1.2
            if "strEdgeAttr" in data:
                nxGraph[3][2]["strEdgeAttr"] = "32"
            if "complexEdgeAttr" in data:
                nxGraph[3][2]["complexEdgeAttr"] = (3, 2, "32")

        return nxGraph

    def test_nx2nk_weight_undirected(self):
        nxG = self.toyNxgraph(directed=False, data=["weight"])
        nkG = nk.nxadapter.nx2nk(nxG, weightAttr="weight", data=False)

        self.assertEqual(nkG.numberOfNodes(), 4)
        self.assertEqual(nkG.numberOfEdges(), 3)

        self.assertEqual(nkG.weight(0, 1), 2)
        self.assertEqual(nkG.weight(0, 2), 4)
        self.assertEqual(nkG.weight(1, 3), -6)

    def test_nx2nk_weight_directed(self):
        nxG = self.toyNxgraph(directed=True, data=["weight"])
        nkG = nk.nxadapter.nx2nk(nxG, weightAttr="weight", data=False)

        self.assertEqual(nkG.numberOfNodes(), 4)
        self.assertEqual(nkG.numberOfEdges(), 7)

        self.assertEqual(nkG.weight(0, 1), 2)
        self.assertEqual(nkG.weight(0, 2), 4)
        self.assertEqual(nkG.weight(1, 3), -6)
        self.assertEqual(nkG.weight(1, 0), 2)
        self.assertEqual(nkG.weight(2, 0), 4)
        self.assertEqual(nkG.weight(3, 1), -6)
        self.assertEqual(nkG.weight(2, 3), -5)

    @unittest.mock.patch("sys.stderr")
    def test_nx2nk_int_undirected(self, mock_stderr):
        nxG = self.toyNxgraph(directed=False, data=["intEdgeAttr", "intNodeAttr"])
        nkG = nk.nxadapter.nx2nk(nxG, data=True)

        edgeAttr = nkG.getEdgeAttribute("intEdgeAttr", int)
        self.assertEqual(edgeAttr[0, 1], 1)
        self.assertEqual(edgeAttr[0, 2], 2)
        self.assertEqual(edgeAttr[1, 3], -3)

        nodeAttr = nkG.getNodeAttribute("intNodeAttr", int)
        self.assertEqual(nodeAttr[0], 4)
        self.assertEqual(nodeAttr[1], 5)
        self.assertEqual(nodeAttr[3], 6)
        self.assertEqual(nodeAttr[2], 7)

        mock_stderr.write.assert_not_called()

    @unittest.mock.patch("sys.stderr")
    def test_nx2nk_int_directed(self, mock_stderr):
        nxG = self.toyNxgraph(directed=True, data=["intEdgeAttr"])
        nkG = nk.nxadapter.nx2nk(nxG, data=True)

        edgeAttr = nkG.getEdgeAttribute("intEdgeAttr", int)
        self.assertEqual(edgeAttr[0, 1], 1)
        self.assertEqual(edgeAttr[0, 2], 2)
        self.assertEqual(edgeAttr[1, 3], -3)
        self.assertEqual(edgeAttr[1, 0], 1)
        self.assertEqual(edgeAttr[2, 0], 2)
        self.assertEqual(edgeAttr[3, 1], -3)
        self.assertEqual(edgeAttr[2, 3], -1)

        mock_stderr.write.assert_not_called()

    @unittest.mock.patch("sys.stderr")
    def test_nx2nk_float_undirected(self, mock_stderr):
        nxG = self.toyNxgraph(directed=False, data=["floatEdgeAttr", "floatNodeAttr"])
        nkG = nk.nxadapter.nx2nk(nxG, data=True)

        edgeAttr = nkG.getEdgeAttribute("floatEdgeAttr", float)
        self.assertEqual(edgeAttr[0, 1], 1.2)
        self.assertEqual(edgeAttr[0, 2], 2.3)
        self.assertEqual(edgeAttr[1, 3], -3.1)

        nodeAttr = nkG.getNodeAttribute("floatNodeAttr", float)
        self.assertEqual(nodeAttr[0], 4.7)
        self.assertEqual(nodeAttr[1], 5.6)
        self.assertEqual(nodeAttr[3], 6.5)
        self.assertEqual(nodeAttr[2], 7.4)

        mock_stderr.write.assert_not_called()

    @unittest.mock.patch("sys.stderr")
    def test_nx2nk_float_directed(self, mock_stderr):
        nxG = self.toyNxgraph(directed=True, data=["floatEdgeAttr"])
        nkG = nk.nxadapter.nx2nk(nxG, data=True)

        edgeAttr = nkG.getEdgeAttribute("floatEdgeAttr", float)
        self.assertEqual(edgeAttr[0, 1], 1.2)
        self.assertEqual(edgeAttr[0, 2], 2.3)
        self.assertEqual(edgeAttr[1, 3], -3.1)
        self.assertEqual(edgeAttr[1, 0], 1.2)
        self.assertEqual(edgeAttr[2, 0], 2.3)
        self.assertEqual(edgeAttr[3, 1], -3.1)
        self.assertEqual(edgeAttr[2, 3], -1.2)

        mock_stderr.write.assert_not_called()

    @unittest.mock.patch("sys.stderr")
    def test_nx2nk_str_undirected(self, mock_stderr):
        nxG = self.toyNxgraph(directed=False, data=["strEdgeAttr", "strNodeAttr"])
        nkG = nk.nxadapter.nx2nk(nxG, data=True)

        edgeAttr = nkG.getEdgeAttribute("strEdgeAttr", str)
        self.assertEqual(edgeAttr[0, 1], "0a")
        self.assertEqual(edgeAttr[0, 2], "03")
        self.assertEqual(edgeAttr[1, 3], "a2")

        nodeAttr = nkG.getNodeAttribute("strNodeAttr", str)
        self.assertEqual(nodeAttr[0], "n0")
        self.assertEqual(nodeAttr[1], "na")
        self.assertEqual(nodeAttr[3], "n2")
        self.assertEqual(nodeAttr[2], "n3")

        mock_stderr.write.assert_not_called()

    @unittest.mock.patch("sys.stderr")
    def test_nx2nk_str_directed(self, mock_stderr):
        nxG = self.toyNxgraph(directed=True, data=["strEdgeAttr"])
        nkG = nk.nxadapter.nx2nk(nxG, data=True)

        edgeAttr = nkG.getEdgeAttribute("strEdgeAttr", str)
        self.assertEqual(edgeAttr[0, 1], "0a")
        self.assertEqual(edgeAttr[0, 2], "03")
        self.assertEqual(edgeAttr[1, 3], "a2")
        self.assertEqual(edgeAttr[1, 0], "0a")
        self.assertEqual(edgeAttr[2, 0], "03")
        self.assertEqual(edgeAttr[3, 1], "a2")
        self.assertEqual(edgeAttr[2, 3], "32")

        mock_stderr.write.assert_not_called()

    def test_nx2nk_complex_undirected(self):
        nxG = self.toyNxgraph(
            directed=False, data=["complexEdgeAttr", "complexNodeAttr"]
        )

        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            nkG = nk.nxadapter.nx2nk(nxG, data=True)
            # Verify some things
            self.assertEqual(len(w), 2)
            self.assertIn(
                "Info: the node attribute complexNodeAttr has been converted to its string representation.",
                str(w[0].message),
            )
            self.assertIn(
                "Info: the edge attribute complexEdgeAttr has been converted to its string representation.",
                str(w[1].message),
            )

        edgeAttr = nkG.getEdgeAttribute("complexEdgeAttr", str)
        self.assertEqual(edgeAttr[0, 1], "(0, 'a', '0a')")
        self.assertEqual(edgeAttr[0, 2], "(0, 3, '03')")
        self.assertEqual(edgeAttr[1, 3], "('a', 2, 'a2')")

        nodeAttr = nkG.getNodeAttribute("complexNodeAttr", str)
        self.assertEqual(nodeAttr[0], "(0, '0')")
        self.assertEqual(nodeAttr[1], "('a', 'a')")
        self.assertEqual(nodeAttr[3], "(2, '2')")
        self.assertEqual(nodeAttr[2], "(3, '3')")

    def test_nx2nk_complex_directed(self):
        nxG = self.toyNxgraph(directed=True, data=["complexEdgeAttr"])

        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            nkG = nk.nxadapter.nx2nk(nxG, data=True)
            # Verify some things
            self.assertEqual(len(w), 1)
            self.assertIn(
                "Info: the edge attribute complexEdgeAttr has been converted to its string representation.",
                str(w[0].message),
            )

        edgeAttr = nkG.getEdgeAttribute("complexEdgeAttr", str)
        self.assertEqual(edgeAttr[0, 1], "(0, 'a', '0a')")
        self.assertEqual(edgeAttr[0, 2], "(0, 3, '03')")
        self.assertEqual(edgeAttr[1, 3], "('a', 2, 'a2')")
        self.assertEqual(edgeAttr[1, 0], "(0, 'a', '0a')")
        self.assertEqual(edgeAttr[2, 0], "(0, 3, '03')")
        self.assertEqual(edgeAttr[3, 1], "('a', 2, 'a2')")
        self.assertEqual(edgeAttr[2, 3], "(3, 2, '32')")

    def test_nx2nk(self):
        nkG = nk.Graph(3)
        nkG.addEdge(0, 1)
        nkG.addEdge(0, 2)

        nxG = nk.nxadapter.nk2nx(nkG)

        self.assertEqual(len(nxG.edges()), 2)
        self.assertEqual(len(nxG.nodes()), 3)

    @unittest.mock.patch("sys.stderr")
    def test_nx2nk_custom_typemap(self, mock_stderr):
        nxG = self.toyNxgraph(directed=False, data=["intEdgeAttr", "intNodeAttr"])
        nkG = nk.nxadapter.nx2nk(nxG, data=True, typeMap={"intEdgeAttr": int})

        edgeAttr = nkG.getEdgeAttribute("intEdgeAttr", int)
        self.assertEqual(edgeAttr[0, 1], 1)
        self.assertEqual(edgeAttr[0, 2], 2)
        self.assertEqual(edgeAttr[1, 3], -3)

        nodeAttr = nkG.getNodeAttribute("intNodeAttr", int)
        self.assertEqual(nodeAttr[0], 4)
        self.assertEqual(nodeAttr[1], 5)
        self.assertEqual(nodeAttr[3], 6)
        self.assertEqual(nodeAttr[2], 7)

        mock_stderr.write.assert_not_called()


if __name__ == "__main__":
    unittest.main()
