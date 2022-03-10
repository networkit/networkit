""" This module provides convenient workflows constructed from NetworKit functions."""

__author__ = "Christian Staudt"

# internal imports
from .graphio import readGraph
from .components import ConnectedComponents
from . import stopwatch 

# external imports
import operator
import logging
import unittest
import os
import csv
import fnmatch

from warnings import warn

def extractLargestComponent(G):
	"""
	extractLargestComponent(G)

	Extract the subgraph of the largest connected component.

	DEPRECATED. This function (and the networkit.workflows module) will be removed in future updates. 

	Parameters
	----------
	G : networkit.Graph
		Input graph.

	Returns
	-------
	networkit.Graph
		Subgraph of largest component, preserving node ids of orignal graph.
	"""

	from warnings import warn
	warn("This function is deprecated, use extractLargestConnectedComponent in the ConnectedComponents module instead.")
	return ConnectedComponents.extractLargestConnectedComponent(G)


def batch(graphDir, match, format, function, outPath, header=None):
	"""
	batch(graphDir, match, format, function, outPath, header=None)

	Read graphs from a directory, apply a function and store result in CSV format.

	DEPRECATED. This function (and the networkit.workflows module) will be removed in future updates. 

	Parameters
	----------
	graphDir : str
		A directory containing graph files.
	match : str
		A pattern that must match the filename so the file is treated as a graph.
	format : networkit.graphio.Format
		Graph file format.
	function : networkit.base.Algorithm
		Any algorithm with a graph as input and a list/tuple of values as output.
	outPath : str
		Path of output CSV file.
	header : str, optional
		CSV file header. Default: None
	"""
	warn("networkit.workflows.batch is deprecated, will be removed in future updates.")
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
					G = readGraph(graphPath)
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
