""" This module provides convenient workflows constructed from NetworKit functions."""

__author__ = "Christian Staudt"

# external imports
import operator
import logging
import unittest
import os
import csv
import fnmatch

import networkit as nk

def extractLargestComponent(G):
	"""
	Extract the subgraph of the largest connected component.

	Parameters:
	-----------
	G : Graph
		Input graph.

	Returns:
	--------
	Graph
		Subgraph of largest component, preserving node ids of orignal graph.
	"""

	from warnings import warn
	warn("This function is deprecated, use extractLargestConnectedComponent in the ConnectedComponents module instead.")
	return nk.components.extractLargestConnectedComponent(G)


def batch(graphDir, match, format, function, outPath, header=None):
	"""
	Read graphs from a directory, apply a function and store result in CSV format.

	Parameters:
	-----------
	graphDir : a directory containing graph files
	match : a pattern that must match the filename so the file is treated as a graph
	format : graph file format
	function : any function from Graph to list/tuple of values
	outPath : path of output CSV file
	header : CSV file header
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
					G = nk.graphio.readGraph(graphPath)
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

if __name__ == '__main__':
    unittest.main()
