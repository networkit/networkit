#!/usr/bin/env python3
import unittest

import networkit as nk

class TestReachability(unittest.TestCase):

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
