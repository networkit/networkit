"""
This module contains general purpose code.
"""
import os
import csv
import fnmatch

import stopwatch
import graphio

# batch processing


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
