#!/usr/bin/env python3
import unittest
import networkit as nk

class TestGraphTools(unittest.TestCase):
	def testLubyAlgorithm(self):
		G = nk.Graph(4, False, False)
		G.addEdge(0, 1)
		G.addEdge(0, 2)
		G.addEdge(1, 2)
		G.addEdge(2, 3)

		luby = nk.independentset.Luby()
		res = luby.run(G)
		count = sum(res)
		# The are several valid outcomes, with either one or two nodes being independent.
		self.assertGreaterEqual(count, 1)

		G.addEdge(0, 3)
		G.addEdge(1, 3)
		res = luby.run(G)
		count = sum(res)
		# Only a single node can be independent since the graph is fully connected.
		self.assertEqual(count, 1)

	def testIndependentSet(self):
		G = nk.Graph(3, False, False)
		G.addEdge(0, 1)
		G.addEdge(0, 2)

		luby = nk.independentset.Luby()
		#Set {0} is independent
		self.assertTrue(luby.isIndependentSet([True,False,False], G))
		#Set {1,2} is independent
		self.assertTrue(luby.isIndependentSet([False,True,True], G))	

if __name__ == "__main__":
	unittest.main()
