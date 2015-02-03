import unittest
import os
#from .. import networkit
#from .. import *
#import ..


class TestExtMETISGraphReader(unittest.TestCase):
	def setUp(self):
		from _NetworKit import Graph
		self.g = Graph(5)
		self.g.addEdge(0,1)
		self.g.addEdge(0,2)
		self.g.addEdge(0,3)
		self.g.addEdge(0,4)

	def test_readAndWrite(self):
#		from networkit.graphio import METISGraphReader
		from _NetworKit import METISGraphReader
		from _NetworKit import METISGraphWriter
		w = METISGraphWriter()
		w.write(self.g, "output/metis_test.graph")
		self.assertTrue(os.path.isfile("output/metis_test.graph"))
		r = METISGraphReader()
		testg = r.read("output/metis_test.graph")
		self.assertEqual(self.g.numberOfNodes(), testg.numberOfNodes())
		self.assertEqual(self.g.numberOfEdges(), testg.numberOfEdges())

