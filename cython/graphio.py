from _NetworKit import (METISGraphReader, FastMETISGraphReader, DotGraphWriter, EdgeListIO, \
						 LineFileReader, SNAPGraphWriter, ClusteringReader, ClusteringWriter)

def readGraph(path, format=None, fast=False):
	"""    Read graph and return a NetworKit::Graph"""
	# TODO: detect file format by looking at the file content
	if format is None:    # look at file extension
		if path.endswith(".graph") or path.endswith(".metis") or path.endswith(".dimacs"):
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
			reader = EdgeListIO('\t', 1)

	with open(path, "r") as file:    # catch a wrong path before it crashes the interpreter
		G = reader.read(path)
		return G

	return None

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
	
	readers = {"metis": METISGraphReader, "edgelist" : EdgeListIO}    
	writers = {"edgelist": EdgeListIO}
	
	reader = readers[fromFormat]()
	writer = writers[toFormat]()
	
	return GraphConverter(reader, writer)