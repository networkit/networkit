""" This module provides convenient workflows constructed from NetworKit functions."""

import operator
import logging
import unittest

from networkit import properties, graph, generators

def extractLargestComponent(G):
	"""
	Extract the subgraph of the largest connected component.

	Parameters
	----------
	G : Graph
		Input graph.
	Returns
	-------
	Graph
		Subgraph of largest component, preserving node ids of orignal graph.
	"""

	cc = properties.ConnectedComponents(G)
	cc.run()
	cSizes = cc.getComponentSizes()
	(largestCompo, size) = max(cSizes.items(), key=operator.itemgetter(1))
	logging.info("extracting component {0} containing {1} nodes".format(largestCompo, size))
	compoNodes = [v for v in G.nodes() if cc.componentOfNode(v) is largestCompo]
	C = graph.Subgraph().fromNodes(G, compoNodes)
	return C



class TestWorkflows(unittest.TestCase):

	def testExtractLargestComponent(self):
		G = generators.DorogovtsevMendesGenerator(100).generate()
		C = extractLargestComponent(G)
		self.assertEqual(properties.size(C), properties.size(G))

if __name__ == '__main__':
    unittest.main()
