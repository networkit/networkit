#!/usr/bin/env python3
import unittest

import networkit as nk

class TestReachability(unittest.TestCase):

	def testAllSimplePaths(self):
		G = nk.Graph(8, directed=True)
		G.addEdge(0,1)
		G.addEdge(2,3)
		G.addEdge(0,5)
		G.addEdge(1,3)
		G.addEdge(4,5)
		G.addEdge(5,6)
		G.addEdge(6,7)
		G.addEdge(3,6)

		asp = nk.reachability.AllSimplePaths(G, 0, 7).run()
		self.assertEqual(asp.numberOfSimplePaths(), 2)
		self.assertIn([0, 1, 3, 6, 7], asp.getAllSimplePaths())
		self.assertIn([0, 5, 6, 7], asp.getAllSimplePaths())

	def testAllSimplePathsForIterator(self):
		G = nk.Graph(4, directed=True)
		G.addEdge(0,1)
		G.addEdge(2,3)
		G.addEdge(0,2)
		G.addEdge(1,3)

		asp = nk.reachability.AllSimplePaths(G, 0, 3).run()
		asp.forAllSimplePaths(lambda x : [0,1,2,3])
		self.assertEqual(asp.numberOfSimplePaths(), 2)

	def testReachableNodes(self):
		for directed in [False, True]:
			for exact in [False, True]:
				g = nk.generators.ErdosRenyiGenerator(100, 0.01, directed).generate()
				rn = nk.reachability.ReachableNodes(g, exact).run()
				for u in g.iterNodes():
					reached = []
					nk.traversal.Traversal.BFSfrom(g, u, lambda v, _: reached.append(v))
					if exact:
						self.assertEqual(rn.numberOfReachableNodes(u), len(reached))
					else:
						self.assertLessEqual(rn.numberOfReachableNodesLB(u), len(reached))
						self.assertGreaterEqual(rn.numberOfReachableNodesUB(u), len(reached))

if __name__ == "__main__":
	unittest.main()
