""" This module provides convenient workflows constructed from NetworKit functions."""

__author__ = "Christian Staudt"

# external imports
import operator
import logging
import unittest
import os
import csv
import fnmatch



from networkit import graph, generators, components

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

	cc = components.ConnectedComponents(G)
	cc.run()
	cSizes = cc.getComponentSizes()
	(largestCompo, size) = max(cSizes.items(), key=operator.itemgetter(1))
	logging.info("extracting component {0} containing {1} nodes".format(largestCompo, size))
	compoNodes = [v for v in G.nodes() if cc.componentOfNode(v) is largestCompo]
	C = G.subgraphFromNodes(compoNodes)
	return C


def batch(graphDir, match, format, function, outPath, header=None):
	"""
	Read graphs from a directory, apply a function and store result in CSV format.
	:param	graphDir	a directory containing graph files
	:param	match		a pattern that must match the filename so the file is treated as a graph
	:param 	format		graph file format
	:param  function	any function from Graph to list/tuple of values
	:param	outPath		path of output CSV file
	:param	header		CSV file header
	"""
	with open(outPath, 'w') as outFile:
		writer = csv.writer(outFile, delimiter='\t')
		if header:
			writer.writerow(header)
		for root, _, filenames in os.walk(graphDir):
			for filename in filenames:
				if fnmatch.fnmatch(filename, match):
					print("processing {0}".format(filename))
					graphPath = os.path.join(root, filename)
					timer = stopwatch.Timer()
					G = graphio.readGraph(graphPath)
					timer.stop()
					result = function(G)
					if type(result) is tuple:
						row = list(result)
					elif type(result) is list:
						row = result
					else:
						row = [result]
					row = [filename, timer.elapsed] + list(row)
					writer.writerow(row)


# TODO: move this to testing module

class TestWorkflows(unittest.TestCase):

	def testExtractLargestComponent(self):
		G = generators.DorogovtsevMendesGenerator(100).generate()
		C = extractLargestComponent(G)
		self.assertEqual(C.size(), G.size())

if __name__ == '__main__':
    unittest.main()
