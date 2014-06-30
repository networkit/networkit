from _NetworKit import (Graph, METISGraphReader, METISGraphWriter, DotGraphWriter, EdgeListWriter, \
						 GMLGraphWriter, LineFileReader, SNAPGraphWriter, DGSWriter, \
						  DGSStreamParser, GraphUpdater, SNAPEdgeListPartitionReader, SNAPGraphReader, EdgeListReader)
from GraphMLIO import GraphMLReader, GraphMLWriter
import os
import logging
import numpy
import scipy.io
from enum import Enum


class AutoNumber(Enum):
	def __new__(cls):
		value = len(cls.__members__) + 1
		obj = object.__new__(cls)
		obj._value_ = value
		return obj


class Format(AutoNumber):
	""" Simple enumeration class to list supported file types """
	SNAP = ()
	EdgeListSpaceZero = ()
	EdgeListSpaceOne = ()
	EdgeListTabZero = ()
	EdgeListTabOne = ()
	METIS = ()
	GraphML = ()
	GML = ()
#	VNA = ()
	EdgeListCommaOne = ()
	GraphViz = ()
#	GDF = ()
	EdgeList = ()
	LFR = ()




# reading

def getReader(fileformat, **kwargs):
	#define your [edgelist] reader here:
	readers =	{
			Format.METIS:			METISGraphReader(),
			Format.GraphML:			GraphMLReader(),
			Format.SNAP:			EdgeListReader('\t',0,'#',False),
			Format.EdgeListCommaOne:	EdgeListReader(',',1,),
			Format.EdgeListSpaceOne:	EdgeListReader(' ',1),
			Format.EdgeListSpaceZero:	EdgeListReader(' ',0),
			Format.EdgeListTabOne:		EdgeListReader('\t',1),
			Format.EdgeListTabZero:		EdgeListReader('\t',0),
			Format.LFR:			EdgeListReader('\t',1)
			}

	try:
		# special case for custom Edge Lists
		if fileformat == Format.EdgeList:
			reader = EdgeListReader(kwargs['separator'],kwargs['firstNode'])
		else:
			reader = readers[fileformat]#(**kwargs)
	except Exception or KeyError:
		raise Exception("unrecognized format/format not supported as input: {0}".format(fileformat))
	return reader


def readGraph(path, fileformat = Format.METIS, **kwargs):
	""" Read graph file in various formats and return a NetworKit::Graph
	    Paramaters: 
		- fileformat: An element of the Format enumeration, default is Format.METIS
		- **kwargs: in case of a custom edge list, provide the defining paramaters as follows:
			"separator=CHAR, firstNode=NODE, commentPrefix=STRING, continuous=BOOL"
			commentPrefix and continuous are optional
	"""
	reader = getReader(fileformat,**kwargs)


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
				raise IOError("{0} is not a valid {1} file: {2}".format(path,fileformat,e))
	return None


def readMat(path):
	""" Reads a Graph from a matlab object file containing an adjacency matrix"""
	matlabObject = scipy.io.loadmat(path)	
	# result is a dictionary of variable names and objects, representing the matlab object
	for (key, value) in matlabObject.items():
		if type(matlabObject[key]) is numpy.ndarray:
			A = matlabObject[key]
			break # found the matrix
			
	(n, n2) = A.shape
	if (n != n2):
		raise Exception("this (%sx%s) matrix is not square".format(n, n2))
	if not ((A.transpose() == A).all()):
		logging.warning("the adjacency matrix is not symmetric")
	G = Graph(n)
	nz = A.nonzero()
	edges = [(u,v) for (u,v) in zip(nz[0], nz[1])]
	for (u,v) in edges:
		G.addEdge(u, v)
	return G


# writing
def getWriter(fileformat):
	writers =	{
			Format.METIS:			METISGraphWriter(),
			Format.GraphML:			GraphMLWriter(),
#			Format.SNAP:			EdgeListWriter('\t',0,'#',False),
			Format.EdgeListCommaOne:	EdgeListWriter(',',1,),
			Format.EdgeListSpaceOne:	EdgeListWriter(' ',1),
			Format.EdgeListSpaceZero:	EdgeListWriter(' ',0),
			Format.EdgeListTabOne:		EdgeListWriter('\t',1),
			Format.EdgeListTabZero:		EdgeListWriter('\t',0),
			Format.GraphViz:		DotGraphWriter(),
			Format.GML:			GMLGraphWriter(),
			Format.LFR:			EdgeListWriter('\t',1)
#			Format.GDF:			GDFGraphWriter(),
#			Format.VNA:			VNAGraphWriter(),
			}
	try:
		# special case for custom Edge Lists
		if fileformat == Format.EdgeList:
			writer = EdgeListIO(kwargs['separator'],kwargs['firstNode'])
		else:
			writer = writers[fileformat]#(**kwargs)
	except KeyError:
		raise Exception("format {0} currently not supported".format(fileformat))
	return writer
def writeGraph(G, path, fileformat = Format.METIS):
	""" Write graph to various output formats. 
		Default format is METIS."""
	writer = getWriter(fileformat)
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

