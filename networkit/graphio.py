# extension imports
from _NetworKit import (METISGraphReader, METISGraphWriter, DotGraphWriter, EdgeListWriter, \
						GMLGraphWriter, LineFileReader, SNAPGraphWriter, DGSWriter, GraphToolBinaryWriter, \
						GraphToolBinaryReader, DGSStreamParser, GraphUpdater, SNAPEdgeListPartitionReader, \
						SNAPGraphReader, EdgeListReader, CoverReader, CoverWriter, EdgeListCoverReader, \
						KONECTGraphReader, GMLGraphReader, MultipleEdgesHandling, ThrillGraphBinaryReader, \
						ThrillGraphBinaryWriter)
from _NetworKit import Graph as __Graph
# local imports
from .GraphMLIO import GraphMLReader, GraphMLWriter
from .GEXFIO import GEXFReader, GEXFWriter
from . import algebraic

# external imports
import os
import logging
try:
	import numpy
except ImportError:
	print("module 'numpy' not available - some functionality will be restricted")
try:
	import scipy.io
except ImportError:
	print("module 'scipy' not available - some functionality will be restricted")
import fnmatch

try:
	from enum import Enum

	class __AutoNumber(Enum):
		def __new__(cls):
			value = len(cls.__members__) + 1
			obj = object.__new__(cls)
			obj._value_ = value
			return obj


	class Format(__AutoNumber):
		""" Simple enumeration class to list supported file types. Currently supported
		file types: SNAP, EdgeListSpaceZero, EdgeListSpaceOne, EdgeListTabZero, EdgeListTabOne,
		METIS, GraphML, GEXF, GML, EdgeListCommaOne, GraphViz, DOT, EdgeList, LFR, KONEC, GraphToolBinary"""
		SNAP = ()
		EdgeListSpaceZero = ()
		EdgeListSpaceOne = ()
		EdgeListTabZero = ()
		EdgeListTabOne = ()
		METIS = ()
		GraphML = ()
		GEXF = ()
		GML = ()
		EdgeListCommaOne = ()
		GraphViz = ()
		DOT = ()
		EdgeList = ()
		LFR = ()
		KONECT = ()
		GraphToolBinary = ()
		MAT = ()
		ThrillBinary = ()

except ImportError:
	print("Update to Python >=3.4 recommended - support for < 3.4 may be discontinued in the future")
	class Format:
		SNAP = "snap"
		EdgeListTabOne = "edgelist-t1"
		EdgeListTabZero = "edgelist-t0"
		EdgeListSpaceOne = "edgelist-s1"
		EdgeListSpaceZero = "edgelist-s0"
		METIS = "metis"
		GraphML = "graphml"
		GEXF = "gexf"
		GML = "gml"
		EdgeListCommaOne = "edgelist-cs1"
		GraphViz = "dot"
		DOT = "dot"
		EdgeList = "edgelist"
		LFR = "edgelist-t1"
		KONECT = "konect"
		GraphToolBinary = "gtbin"
		MAT = "mat"
		ThrillBinary = "thrillbinary"





# reading

def getReader(fileformat, **kwargs):
	#define your [edgelist] reader here:
	readers =	{
			Format.METIS:			METISGraphReader(),
			Format.GraphML:			GraphMLReader(),
			Format.GEXF:			GEXFReader(),
			Format.SNAP:			EdgeListReader('\t',0,'#',False),
			Format.EdgeListCommaOne:	EdgeListReader(',',1,),
			Format.EdgeListSpaceOne:	EdgeListReader(' ',1),
			Format.EdgeListSpaceZero:	EdgeListReader(' ',0),
			Format.EdgeListTabOne:		EdgeListReader('\t',1),
			Format.EdgeListTabZero:		EdgeListReader('\t',0),
			Format.LFR:			EdgeListReader('\t',1),
			Format.KONECT:			KONECTGraphReader(),
			Format.GML:			GMLGraphReader(),
			Format.GraphToolBinary:		GraphToolBinaryReader(),
			Format.MAT:			MatReader(),
			Format.ThrillBinary:		ThrillGraphBinaryReader(),
			}

	try:
		# special case for custom Edge Lists
		if fileformat == Format.EdgeList:
			if kwargs["continuous"] == False:
				kwargs["firstNode"] = 0
			reader = EdgeListReader(**kwargs)
		else:
			reader = readers[fileformat]#(**kwargs)
	except Exception or KeyError:
		raise Exception("unrecognized format/format not supported as input: {0}".format(fileformat))
	return reader


def readGraph(path, fileformat, **kwargs):
	""" Read graph file in various formats and return a NetworKit::Graph
	    Parameters:
		- fileformat: An element of the Format enumeration. Currently supported file types:
		SNAP, EdgeListSpaceZero, EdgeListSpaceOne, EdgeListTabZero, EdgeListTabOne, METIS,
		GraphML, GEXF, GML, EdgeListCommaOne, GraphViz, DOT, EdgeList, LFR, KONECT, GraphToolBinary, ThrillBinary
		- **kwargs: in case of a custom edge list, pass the genereic Fromat.EdgeList accompanied by
			the defining paramaters as follows:
			"separator=CHAR, firstNode=NODE, commentPrefix=STRING, continuous=BOOL, directed=BOOL"
			commentPrefix='#', continuous=True and directed=False are optional because of their default values;
			firstNode is not needed when continuous=True.
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
				G.setName(os.path.basename(path).split(".")[0])	# set name of graph to name of file
				return G
			except Exception as e:
				raise IOError("{0} is not a valid {1} file: {2}".format(path,fileformat,e))
	return None

def readGraphs(dirPath, pattern, fileformat, some=None, exclude=None, **kwargs):
	"""
	Read all graph files contained in a directory whose filename contains the pattern, return a dictionary of name to Graph object.
    Parameters:
	- pattern: unix-style string pattern
	- fileformat: An element of the Format enumeration
	- some: restrict number of graphs to be read
	- **kwargs: in case of a custom edge list, provide the defining paramaters as follows:
		"separator=CHAR, firstNode=NODE, commentPrefix=STRING, continuous=BOOL"
		commentPrefix and continuous are optional
	"""
	graphs = {}
	for root, dirs, files in os.walk(dirPath):
		for file in files:
			if fnmatch.fnmatch(file, pattern):
				if (exclude is None) or (not fnmatch.fnmatch(file, exclude)):
					G = readGraph(os.path.join(root, file), fileformat, **kwargs)
					graphs[G.getName()] = G
					if some:
						if len(graphs) == some:
							return graphs
	return graphs


class MatReader:
	def __init__(self, key = "G"):
		self.key = key

	def read(self, path):
		return readMat(path, self.key)

def readMat(path, key="G"):
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
#	if not numpy.array_equal(A, A.transpose): # FIXME this is slow and doesn't work as expected, seems to be False for valid inputs
#		logging.warning("the adjacency matrix is not symmetric")
	G = __Graph(n)
	nz = A.nonzero()
	for (u,v) in zip(nz[0], nz[1]):
		if not G.hasEdge(u, v):
			G.addEdge(u, v)
	return G

class MatWriter:
	def __init__(self):
		self.key = key

	def write(self, G, path, key="G"):
		writeMat(path, key)

def writeMat(G, path, key="G"):
	""" Writes a NetworKit::Graph to a Matlab object file.
		Parameters:
		- G: The graph
		- path: Path of the matlab file
		- key: Dictionary Key
	"""
	matrix = algebraic.adjacencyMatrix(G, matrixType='sparse')
	scipy.io.savemat(path, {key : matrix})


# writing
def getWriter(fileformat, **kwargs):
	writers =	{
			Format.METIS:			METISGraphWriter(),
			Format.GraphML:			GraphMLWriter(),
			Format.GEXF:			GEXFWriter(),
#			Format.SNAP:			EdgeListWriter('\t',0,'#',False),
			Format.EdgeListCommaOne:	EdgeListWriter(',',1,),
			Format.EdgeListSpaceOne:	EdgeListWriter(' ',1),
			Format.EdgeListSpaceZero:	EdgeListWriter(' ',0),
			Format.EdgeListTabOne:		EdgeListWriter('\t',1),
			Format.EdgeListTabZero:		EdgeListWriter('\t',0),
			Format.GraphViz:		DotGraphWriter(),
			Format.DOT:			DotGraphWriter(),
			Format.GML:			GMLGraphWriter(),
			Format.LFR:			EdgeListWriter('\t',1),
			Format.GraphToolBinary:		GraphToolBinaryWriter()
			}
	try:
		# special case for custom Edge Lists
		if fileformat == Format.EdgeList:
			writer = EdgeListWriter(**kwargs)
		else:
			writer = writers[fileformat]#(**kwargs)
	except KeyError:
		raise Exception("format {0} currently not supported".format(fileformat))
	return writer

def writeGraph(G, path, fileformat, **kwargs):
	""" Write graph to various output formats.

	Paramaters:
	- G:			a graph
	- path: 		output path
	- fileformat: 	an element of the Format enumeration

	"""

	dirname = os.path.dirname(os.path.realpath(path))
	# the given file path does not exist yet
	if not os.path.isfile(path):
		# check write permissions on the directory
		if not os.access(dirname, os.W_OK):
			# we may not write on this directory, raise Error
			raise IOError("No permission to write")
		# else everthing is alright
	else:
		# the given path points to a file
		if not os.access(path, os.W_OK):
			raise IOError("No permission to write")
		else:
			logging.warning("overriding given file")
	writer = getWriter(fileformat, **kwargs)
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
