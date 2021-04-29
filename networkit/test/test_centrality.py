#!/usr/bin/env python3
import unittest
import os
import networkit as nk


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

if __name__ == "__main__":
	unittest.main()
