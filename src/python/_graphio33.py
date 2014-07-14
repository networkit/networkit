from _NetworKit import (Graph, METISGraphReader, METISGraphWriter, DotGraphWriter, EdgeListWriter, \
						 GMLGraphWriter, LineFileReader, SNAPGraphWriter, DGSWriter, \
						  DGSStreamParser, GraphUpdater, SNAPEdgeListPartitionReader, SNAPGraphReader, EdgeListReader, CoverReader)
from GraphMLIO import GraphMLReader, GraphMLWriter
import os
import logging
import numpy
import scipy.io

class Format:
	METIS = "metis"
	GraphML = "graphml"
	EdgeListTabOne = "edgelist-t1"
	EdgeListTabZero = "edgelist-t0"
	EdgeListSpaceOne = "edgelist-s1"
	EdgeListSpaceZero = "edgelist-s0"
	EdgeListCommaOne = "edgelist-cs1"
	SNAP = "snap"
	GraphViz = "dot"
	DOT = "dot"
	EdgeList = "edgelist"
	LFR = "edgelist-t1"
	GML = "gml"


# reading
def getReader(format, **kwargs):
	readers =  {
			"metis": METISGraphReader(),
			"edgelist-t1": EdgeListReader('\t', 1),
			"edgelist-t0": EdgeListReader('\t', 0),
			"edgelist-s1": EdgeListReader(' ', 1),
			"edgelist-s0": EdgeListReader(' ', 0),
			"edgelist-cs1": EdgeListReader(',',1),
			"graphml": GraphMLReader(),
			"snap": EdgeListReader('\t',0,'#',False),
		}
	try:
		# special case for custom Edge Lists
		if format == "edgelist":
			reader = EdgeListReader(kwargs['separator'],kwargs['firstNode'])
		else:
			reader = readers[format]#(**kwargs)
	except Exception or KeyError:
		raise Exception("unrecognized format/format not supported as input: {0}".format(format))
	return reader


def readGraph(path, format="metis", **kwargs):
	"""    Read graph file in various formats and return a NetworKit::Graph
		Default format is METIS"""
	reader = getReader(format, **kwargs)

	if ("~" in path):
		path = os.path.expanduser(path)
		print("path expanded to: {0}".format(path))
	if not os.path.isfile(path):
		raise IOError("{0} is not a file".format(path))
	else:
		with open(path, "r") as file:    # catch a wrong path before it crashes the interpreter
			try:
				G = reader.read(path)
				return G
			except Exception as e:
				raise IOError("{0} is not a valid {1} file: {2}".format(path,format,e))
	return None


def readMat(path, key="A"):
	""" Reads a Graph from a matlab object file containing an adjacency matrix and returns a NetworKit::Graph
		Parameters:
		- key: The key of the adjacency matrix in the matlab object file (default: A)"""
	matlabObject = scipy.io.loadmat(path)
	# result is a dictionary of variable names and objects, representing the matlab object
	if key in matlabObject:
		A = matlabObject[key]
	else:
		raise Exception("Key {0} not found in the matlab object file".format(key))
	(n, n2) = A.shape
	if n != n2:
		raise Exception("this ({0}x{1}) matrix is not square".format(n, n2))
	if not numpy.array_equal(A, A.transpose):
		logging.warning("the adjacency matrix is not symmetric")
	G = Graph(n)
	nz = A.nonzero()
	for (u,v) in zip(nz[0], nz[1]):
		G.addEdge(u, v)
	return G

# writing
def getWriter(format, **kwargs):
	writers = {
			"metis" : METISGraphWriter(),
			"gexf": None,
#			"vna": VNAGraphWriter(),
			"dot": DotGraphWriter(),
			"graphviz": DotGraphWriter(),
			"gml": GMLGraphWriter(),
			"edgelist-t1" : EdgeListWriter('\t', 1),
			"edgelist-t0": EdgeListWriter('\t', 0),
			"edgelist-s1": EdgeListWriter(' ', 1),
			"edgelist-s0": EdgeListWriter(' ', 1),
			"edgelist-cs1": EdgeListWriter(',',1),
			"graphml": GraphMLWriter()
		}
	try:
		# special case for custom Edge Lists
		if format == "edgelist":
			writer = EdgeListWriter(kwargs['separator'],kwargs['firstNode'])
		else:
			writer = writers[format]#(**kwargs)
	except Exception or KeyError:
		raise Exception("format {0} currently not supported".format(format))
	return writer

def writeGraph(G, path, format="metis", **kwargs):
	""" Write graph to various output formats """
	writer = getWriter(format, **kwargs)
	writer.write(G, path)
	logging.info("wrote graph {0} to file {1}".format(G, path))


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

	reader = getReader(fromFormat)
	writer = getWriter(toFormat)

	return GraphConverter(reader, writer)


def convertGraph(fromFormat, toFormat, fromPath, toPath=None):
	converter = getConverter(fromFormat, toFormat)
	if toPath is None:
		toPath = "{0}.{1}.graph".format(fromPath.split(".")[0], toFormat)
	converter.convert(fromPath, toPath)
	print("converted {0} to {1}".format(fromPath, toPath))



# dynamic

def readStream(path, mapped=True, baseIndex=0):
	"""
		Read a graph event stream from a file.
	"""
	return DGSStreamParser(path, mapped, baseIndex).getStream()

def writeStream(stream, path):
	"""
		Write a graph event stream to a file.
	"""
	DGSWriter().write(stream, path)


def graphFromStreamFile(path, mapped=True, baseIndex=0):
	stream = readStream(path, mapped, baseIndex)
	G = Graph()
	gu = GraphUpdater(G)
	gu.update(stream)
	return G
