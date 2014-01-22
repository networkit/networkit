from _NetworKit import (METISGraphReader, FastMETISGraphReader, METISGraphWriter, DotGraphWriter, EdgeListIO, \
						 LineFileReader, SNAPGraphWriter, ClusteringReader, ClusteringWriter)

import os
import logging

def readGraph(path, format=None, fast=False, separator='\t', firstNode=1):
	"""    Read graph and return a NetworKit::Graph"""
	# TODO: detect file format by looking at the file content
	if format is None:    # look at file extension
		if path.endswith(".graph") or path.endswith(".metis") or path.endswith(".dimacs") or path.endswith(".metis.graph"):
			if fast:
				print("using experimental fast reader")
				reader = FastMETISGraphReader()
			else:
				reader = METISGraphReader()
		else:
			raise Exception("unknown graph file format")
	else:
		if (format is "metis"):
			if fast:
				print("using experimental fast reader")
				reader = FastMETISGraphReader()
			else:
				reader = METISGraphReader()
		elif (format is "edgelist"):
			reader = EdgeListIO(separator, firstNode)

	if ("~" in path):
		path = os.path.expanduser(path)
		print("path expanded to: {0}".format(path))
	if not os.path.isfile(path):
		raise IOError("{0} is not a file".format(path))
	else:
		with open(path, "r") as file:    # catch a wrong path before it crashes the interpreter
			G = reader.read(path)
			return G

	return None


def writeGraph(G, path, format="metis"):
	if format is "metis":
		writer = METISGraphWriter()
		writer.write(G, path)
		logging.info("wrote graph {0} to file {1}".format(G, path))
	else:
		raise Exception("format {0} currently not supported".format(format))		

class GraphConverter:
	
	def __init__(self, reader, writer):
		self.reader = reader
		self.writer = writer
		
	def convert(self, inPath, outPath):
		G = self.reader.read(inPath)
		self.writer.write(G, outPath)
		
	def __str__(self):
		return "GraphConverter: {0} => {0}".format(self.reader, self.writer)

def getConverter(fromFormat, toFormat):
	
	readers = {"metis": METISGraphReader(),
				"edgelist-t1" : EdgeListIO('\t', 1),
				 "edgelist-t0": EdgeListIO('\t', 0),
				 "edgelist-s1": EdgeListIO(' ', 1), 
				 "edgelist-s0": EdgeListIO(' ', 1)}    
	writers = {"metis": METISGraphWriter(),
				"edgelist-t1" : EdgeListIO('\t', 1),
				 "edgelist-t0": EdgeListIO('\t', 0),
				 "edgelist-s1": EdgeListIO(' ', 1), 
				 "edgelist-s0": EdgeListIO(' ', 1)}   
	
	reader = readers[fromFormat]
	writer = writers[toFormat]
	
	return GraphConverter(reader, writer)


def convertGraph(fromFormat, toFormat, fromPath, toPath=None):
	converter = getConverter(fromFormat, toFormat)
	if toPath is None:
		toPath = "{0}.{1}.graph".format(fromPath.split(".")[0], toFormat)
	converter.convert(fromPath, toPath)
	print("converted {0} to {1}".format(fromPath, toPath))
