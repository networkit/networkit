# distutils: language=c++

from libc.stdint cimport uint64_t
from libc.stdint cimport uint8_t

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.map cimport map
from libcpp.unordered_map cimport unordered_map

import os
import logging
import numpy
import scipy.io
import fnmatch
from enum import Enum

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node

from .dynamics import DGSWriter, DGSStreamParser
from .graph cimport _Graph, Graph
from .graph import Graph as __Graph
from .structures cimport _Cover, Cover, _Partition, Partition
from .GraphMLIO import GraphMLReader, GraphMLWriter
from .GEXFIO import GEXFReader, GEXFWriter
from . import algebraic
from .support import MissingDependencyError

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Graph move( _Graph t ) nogil # specialized declaration as general declaration disables template argument deduction and doesn't work
	_Partition move( _Partition t) nogil

def stdstring(pystring):
	""" 
	stdstring(pystring)

	Convert a Python string to a bytes object which is automatically coerced to std::string

	Parameters
	----------
	pystring : str
		Input python string.

	Returns
	-------
	stdstring
		Python bytes string.
	"""
	pybytes = pystring.encode("utf-8")
	return pybytes

cdef extern from "<networkit/io/GraphReader.hpp>":

	cdef cppclass _GraphReader "NetworKit::GraphReader":
		_GraphReader() nogil except +
		_Graph read(string path) nogil except +

cdef class GraphReader:
	""" Abstract base class for graph readers"""

	cdef _GraphReader* _this

	def __init__(self, *args, **kwargs):
		if type(self) == GraphReader:
			raise RuntimeError("Error, you may not use GraphReader directly, use a sub-class instead")

	def __cinit__(self, *args, **kwargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

	def read(self, path):
		"""
		read(path)

		Read graph given by path.

		Parameters
		----------
		path : str
			Path string.

		Returns
		-------
		networkit.Graph
			The resulting graph.
		"""
		cdef string cpath = stdstring(path)
		cdef _Graph result

		with nogil:
			result = move(self._this.read(cpath)) # extra move in order to avoid copying the internal variable that is used by Cython
		return Graph(0).setThis(result)

cdef extern from "<networkit/io/GraphReader.hpp>" namespace "NetworKit::GraphReader":

	cdef enum _MultipleEdgesHandling "NetworKit::GraphReader::MultipleEdgesHandling":
		DISCARD_EDGES,
		SUM_WEIGHTS_UP,
		KEEP_MINIMUM_WEIGHT

class MultipleEdgesHandling:
	DiscardEdges = DISCARD_EDGES
	SumWeightsUp = SUM_WEIGHTS_UP
	KeepMinimumWeight = KEEP_MINIMUM_WEIGHT

cdef extern from "<networkit/io/GraphWriter.hpp>":

	cdef cppclass _GraphWriter "NetworKit::GraphWriter":
		_GraphWriter() nogil except +
		void write(_Graph G, string path) nogil except +


cdef class GraphWriter:
	"""
	Abstract base class for graph writers
	"""

	cdef _GraphWriter *_this

	def __init__(self, *args, **kwargs):
		if type(self) == GraphWriter:
			raise RuntimeError("Error, you may not use GraphWriter directly, use a sub-class instead")

	def __cinit__(self, *args, **kwargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

	def write(self, Graph G not None, path):
		"""
		Write the graph to a file.

		Parameters
		----------
		G : networkit.Graph
			The graph to write.
		paths : str
			The output path.
		"""
		assert path != None
		cdef string c_path = stdstring(path)
		with nogil:
			self._this.write(G._this, c_path)
		return self

cdef extern from "<networkit/io/METISGraphReader.hpp>":

	cdef cppclass _METISGraphReader "NetworKit::METISGraphReader" (_GraphReader):
		_METISGraphReader() nogil except +

cdef class METISGraphReader(GraphReader):
	"""
	METISGraphReader()
	
	Reads the METIS adjacency file format [1]. If the Fast reader fails,
	use readGraph(path, graphio.formats.metis) as an alternative.
	[1]: http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
	"""
	def __cinit__(self):
		self._this = new _METISGraphReader()

cdef extern from "<networkit/io/NetworkitBinaryReader.hpp>":
	cdef cppclass _NetworkitBinaryReader "NetworKit::NetworkitBinaryReader" (_GraphReader):
		_NetworkitBinaryReader() except +
		_Graph readFromBuffer(vector[uint8_t] state) except +

cdef class NetworkitBinaryReader(GraphReader):
	"""
	NetworkitBinaryReader()

	Reads a graph written in the custom Networkit format. Further information can be found here: 
	https://github.com/networkit/networkit/blob/master/networkit/cpp/io/NetworkitBinaryGraph.md
	"""

	def __cinit__(self):
		self._this = new _NetworkitBinaryReader()
	
	def readFromBuffer(self, state):
		"""
		readFromBuffer(state)

		Read graph based on input buffer.

		Parameters
		----------
		buffer : list(int)
			Input data buffer.
		"""
		cdef _Graph result
		result = move((<_NetworkitBinaryReader*>(self._this)).readFromBuffer(state))
		return Graph(0).setThis(result)

cdef extern from "<networkit/io/NetworkitBinaryWriter.hpp>":
	cdef cppclass _NetworkitBinaryWriter "NetworKit::NetworkitBinaryWriter" (_GraphWriter):
		_NetworkitBinaryWriter() except +
		vector[uint8_t] writeToBuffer(_Graph G) except +

cdef class NetworkitBinaryWriter(GraphWriter):
	"""
	NetworkitBinaryWriter()

	Writes a graph written in the custom Networkit format. Further information can be found here:
	https://github.com/networkit/networkit/blob/master/networkit/cpp/io/NetworkitBinaryGraph.md
	"""
	def __cinit__(self):
		self._this = new _NetworkitBinaryWriter()
	
	def writeToBuffer(self, Graph G not None):
		"""
		writeToBuffer(state)

		Write graph to data buffer.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		return (<_NetworkitBinaryWriter*>(self._this)).writeToBuffer(G._this)

cdef extern from "<networkit/io/GraphToolBinaryReader.hpp>":

	cdef cppclass _GraphToolBinaryReader "NetworKit::GraphToolBinaryReader" (_GraphReader):
		_GraphToolBinaryReader() except +

cdef class GraphToolBinaryReader(GraphReader):
	""" 
	GraphToolsBinaryReader()

	Reads the binary file format defined by graph-tool: http://graph-tool.skewed.de/static/doc/gt_format.html
	"""
	def __cinit__(self):
		self._this = new _GraphToolBinaryReader()

cdef extern from "<networkit/io/ThrillGraphBinaryReader.hpp>":

	cdef cppclass _ThrillGraphBinaryReader "NetworKit::ThrillGraphBinaryReader" (_GraphReader):
		_ThrillGraphBinaryReader(count n) except +
		_Graph read(vector[string] paths) nogil except +

cdef class ThrillGraphBinaryReader(GraphReader):
	"""
	ThrillGraphBinaryReader(n=0)

	Reads a graph format consisting of a serialized DIA of vector<uint32_t> from thrill.
	When the number of nodes is given, reading the graph is more efficient.
	Otherwise nodes are added to the graph as they are encountered.
	Edges must be present only in one direction.

	Parameters
	----------
	n : int, optional
		The number of nodes. Default: 0
	"""
	def __cinit__(self, count n = 0):
		self._this = new _ThrillGraphBinaryReader(n)

	def read(self, paths):
		"""
		read(paths)

		Read the graph from one or multiple files.

		Parameters
		----------
		paths : str or list(str)
			The input path(s).

		Returns
		-------
		networkit.Graph
			The resulting graph.
		"""
		cdef vector[string] c_paths

		if isinstance(paths, str):
			c_paths.push_back(stdstring(paths))
		else:
			c_paths.reserve(len(paths))

			for p in paths:
				c_paths.push_back(stdstring(p))

		cdef _Graph result

		with nogil:
			result = move((<_ThrillGraphBinaryReader*>(self._this)).read(c_paths)) # extra move in order to avoid copying the internal variable that is used by Cython

		return Graph(0).setThis(result)

cdef extern from "<networkit/io/ThrillGraphBinaryWriter.hpp>":

	cdef cppclass _ThrillGraphBinaryWriter "NetworKit::ThrillGraphBinaryWriter" (_GraphWriter):
		_ThrillGraphBinaryWriter() except +

cdef class ThrillGraphBinaryWriter(GraphWriter):
	"""
	ThrillGraphBinaryWriter()

	Writes a graph format consisting of a serialized DIA of vector<uint32_t> from Thrill.
	Edges are written only in one direction.
	"""

	def __cinit__(self):
		self._this = new _ThrillGraphBinaryWriter()

cdef extern from "<networkit/io/EdgeListReader.hpp>":

	cdef cppclass _EdgeListReader "NetworKit::EdgeListReader"(_GraphReader):
		_EdgeListReader() except +
		_EdgeListReader(char separator, node firstNode, string commentPrefix, bool_t continuous, bool_t directed)
		map[string,node] getNodeMap() except +

cdef class EdgeListReader(GraphReader):
	""" 
	EdgeListReader(self, separator, firstNode, commentPrefix="#", continuous=True, directed=False)

	Reads a graph from various text-based edge list formats.

	A line has to contain two or three entries separated with the separator symbol (one ASCII character).
	If at least one line contains three entries, the generated graph will be weighted and
	each line with only two fields will be interpreted as weight 1.0.

	A file may contain the same edge multiple times; then, the weight of the first
	occurrence is used.

	Undirected graphs need to include an edge only in one direction, i.e. edge {u, v} may
	be represented by (u, v) or (v, u) or both (again, only the first occurrence is used).

	If the input file contains non-continuous node ids with large gaps or non-integer node labels,
	set the parameter continuous to False. Then, gaps are automatically removed and node ids are
	reassigned to [0, n) where n is the number of nodes in the graph. The mapping will be arbitrary
	and can be accessed using getNodeMap().

	To shift continuous integer node labels which are not zero-indexed, set firstNode to
	the smallest id used in the file. firstNode will be ignored in the non-continuous case.

	The file may also include line comments which start with the commentPrefix.

	Parameters
	----------
	separator : str
		The separator character. Must have length of exactly one.
	firstNode : int
		The id of the first node, this value will be subtracted from all node ids.
	commentPrefix : str, optional
		Lines starting with this prefix will be ignored. Default: `#`
	continuous : bool, optional
		File uses continuous node ids. Default: True 
	directed : bool, optional
		Treat input file as a directed graph. Default: False
	"""
	def __cinit__(self, separator, firstNode, commentPrefix="#", continuous=True, directed=False):
		if len(separator) != 1 or ord(separator[0]) > 255:
			raise RuntimeError("separator has to be exactly one ascii character");

		self._this = new _EdgeListReader(stdstring(separator)[0], firstNode, stdstring(commentPrefix), continuous, directed)

	def getNodeMap(self):
		""" 
		getNodeMap()
		
		Returns mapping of non-continuous files.

		The mapping is returned as dict(str, int) projecting the original
		labels (as Strings) to the reassigned node ids.

		Returns
		-------
		dict(str,int)
			Mapping from labels to node ids.
		"""
		cdef map[string,node] cResult = (<_EdgeListReader*>(self._this)).getNodeMap()
		result = dict()
		for elem in cResult:
			result[(elem.first).decode("utf-8")] = elem.second
		return result

cdef extern from "<networkit/io/KONECTGraphReader.hpp>":

	cdef cppclass _KONECTGraphReader "NetworKit::KONECTGraphReader"(_GraphReader):
		_KONECTGraphReader() except +
		_KONECTGraphReader(bool_t remapNodes, _MultipleEdgesHandling handlingmethod)

cdef class KONECTGraphReader(GraphReader):
	""" 
	KONECTGraphReader(remapNodes = False, handlingmethod = networkit.graphio.MultipleEdgesHandling.DiscardEdges)

	Reader for the KONECT graph format, which is described in detail on the KONECT website: http://konect.uni-koblenz.de/downloads/konect-handbook.pdf

	Parameter :code:`handlingmethod` can be one of the following:

	- networkit.graphio.MultipleEdgesHandling.DiscardEdges
	- networkit.graphio.MultipleEdgesHandling.SumWeightsUp
	- networkit.graphio.MultipleEdgesHandling.KeepMinimumWeight

	Parameters
	----------
	remapNodes : bool, optional
		Indicates whether nodes are remapped. Default: False
	handlingmethod : networkit.graphio.MultipleEdgesHandling, optional
		Sets method of handling multiple edges. Default: networkit.graphio.MultipleEdgesHandling.DiscardEdges
	"""
	def __cinit__(self, remapNodes = False, handlingmethod = MultipleEdgesHandling.DiscardEdges):
		self._this = new _KONECTGraphReader(remapNodes, handlingmethod)

cdef extern from "<networkit/io/GMLGraphReader.hpp>":

	cdef cppclass _GMLGraphReader "NetworKit::GMLGraphReader"(_GraphReader):
		_GMLGraphReader() except +

cdef class GMLGraphReader(GraphReader):
	""" 
	GMLGraphReader()

	Reader for the GML graph format, which is documented here: http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf
 	"""
	def __cinit__(self):
		self._this = new _GMLGraphReader()

cdef extern from "<networkit/io/METISGraphWriter.hpp>":

	cdef cppclass _METISGraphWriter "NetworKit::METISGraphWriter" (_GraphWriter):
		_METISGraphWriter() except +


cdef class METISGraphWriter(GraphWriter):
	""" 
	METISGraphWriter()
	
	Writes graphs in the METIS format.
	"""

	def __cinit__(self):
		self._this = new _METISGraphWriter()
cdef extern from "<networkit/io/GraphToolBinaryWriter.hpp>":

	cdef cppclass _GraphToolBinaryWriter "NetworKit::GraphToolBinaryWriter" (_GraphWriter):
		_GraphToolBinaryWriter() except +

cdef class GraphToolBinaryWriter(GraphWriter):
	"""
	GraphToolBinaryWriter()

	Reads the binary file format defined by graph-tool: http://graph-tool.skewed.de/static/doc/gt_format.html
	"""
	def __cinit__(self):
		self._this = new _GraphToolBinaryWriter()

cdef extern from "<networkit/io/DotGraphWriter.hpp>":

	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter" (_GraphWriter):
		_DotGraphWriter() except +

cdef class DotGraphWriter(GraphWriter):
	""" 
	DotGraphWriter()
	
	Writes graphs in the .dot/GraphViz format.
	"""
	def __cinit__(self):
		self._this = new _DotGraphWriter()

cdef extern from "<networkit/io/GMLGraphWriter.hpp>":

	cdef cppclass _GMLGraphWriter "NetworKit::GMLGraphWriter" (_GraphWriter):
		_GMLGraphWriter() except +


cdef class GMLGraphWriter(GraphWriter):
	""" 
	GMLGraphWriter()

	Writes a graph and its coordinates as a GML file: http://svn.bigcat.unimaas.nl/pvplugins/GML/trunk/docs/gml-technical-report.pdf
	"""

	def __cinit__(self):
		self._this = new _GMLGraphWriter()

cdef extern from "<networkit/io/EdgeListWriter.hpp>":

	cdef cppclass _EdgeListWriter "NetworKit::EdgeListWriter" (_GraphWriter):
		_EdgeListWriter() except +
		_EdgeListWriter(char separator, node firstNode, bool_t bothDirections) except +

cdef class EdgeListWriter(GraphWriter):
	""" 
	EdgeListWriter(separator, firstNode, bothDirections = False)

	Writes graphs in various edge list formats.

	Parameters
	----------
	separator : str
		The separator character.
	firstNode : int
		The id of the first node, this value will be added to all node ids
	bothDirections : bool, optional
		If undirected edges shall be written in both directions, i.e., as symmetric directed graph. Default: False
	"""

	def __cinit__(self, separator, firstNode, bool_t bothDirections = False):
		cdef char sep = stdstring(separator)[0]
		self._this = new _EdgeListWriter(sep, firstNode, bothDirections)

cdef extern from "<networkit/io/LineFileReader.hpp>":

	cdef cppclass _LineFileReader "NetworKit::LineFileReader":
		_LineFileReader() except +
		vector[string] read(string path)


cdef class LineFileReader:
	"""
	LineFileReader()
	
	Reads a file and puts each line in a list of strings.
	"""
	cdef _LineFileReader _this

	def read(self, path):
		"""
		read(path)
		
		Reads a file and returns list of strings.
		
		Parameters
		----------
		path : str
			The input path.

		Returns
		-------
		list(str)
			List of strings, each string representing one line of an input file.
		"""
		return self._this.read(stdstring(path))


cdef extern from "<networkit/io/SNAPGraphWriter.hpp>":
	cdef cppclass _SNAPGraphWriter "NetworKit::SNAPGraphWriter" (_GraphWriter):
		_SNAPGraphWriter() except +

cdef class SNAPGraphWriter(GraphWriter):
	"""
	SNAPGraphWriter()
	
	Writes graphs in a format suitable for the Georgia Tech SNAP software: http://snap-graph.sourceforge.net/
	"""

	def __cinit__(self):
		self._this = new _SNAPGraphWriter()

cdef extern from "<networkit/io/SNAPGraphReader.hpp>":

	cdef cppclass _SNAPGraphReader "NetworKit::SNAPGraphReader"(_GraphReader):
		_SNAPGraphReader() except +
		_SNAPGraphReader(bool_t directed, bool_t remapNodes, count nodeCount)

cdef class SNAPGraphReader(GraphReader):
	"""
	SNAPGraphReader(directed = False, remapNodes = True, nodeCount = 0)

	Reads a graph from the SNAP graph data collection: http://snap.stanford.edu/data/index.html

	Parameters
	----------
	directed : bool, optional
		Indicates whether input represents a directed graph. Default: False
	remapNodes : bool, optional
		Indicates whether nodes should be remapped. Default: True
	nodeCount : int, optional
		Indicate the first node id. Default: 0 
	"""
	def __cinit__(self, directed = False, remapNodes = True, nodeCount = 0):
		self._this = new _SNAPGraphReader(directed, remapNodes, nodeCount)



cdef extern from "<networkit/io/PartitionReader.hpp>":

	cdef cppclass _PartitionReader "NetworKit::PartitionReader":
		_PartitionReader() except +
		_Partition read(string path) except +


cdef class PartitionReader:
	""" 
	PartitionReader()

	Reads a partition from a file.
	File format: line i contains subset id of element i.
	 """
	cdef _PartitionReader _this

	def read(self, path):
		"""
		read(path)
		
		Reads a partition from a file.
		
		Parameters
		----------
		path : str
			The input path.

		Returns
		-------
		networkit.Partition
			The resulting partition.
		"""
		return Partition().setThis(self._this.read(stdstring(path)))


cdef extern from "<networkit/io/PartitionWriter.hpp>":

	cdef cppclass _PartitionWriter "NetworKit::PartitionWriter":
		_PartitionWriter() except +
		void write(_Partition, string path) nogil except +


cdef class PartitionWriter:
	"""
	PartitionWriter()
	
	Writes a partition to a file.
	File format: line i contains subset id of element i.
	"""
	cdef _PartitionWriter _this

	def write(self, Partition zeta, path):
		"""
		write(zeta, path)
		
		Writes a partition to a file.
		File format: line i contains subset id of element i.

		Parameters
		----------
		zeta : networkit.Partition
			The input partition.
		path : str
			The output path.
		"""
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(zeta._this, cpath)

cdef extern from "<networkit/io/BinaryPartitionReader.hpp>":

	cdef cppclass _BinaryPartitionReader "NetworKit::BinaryPartitionReader":
		_BinaryPartitionReader() except +
		_BinaryPartitionReader(uint8_t width) except +
		_Partition read(string path) except +


cdef class BinaryPartitionReader:
	"""
	BinaryPartitionReader(width=4)

	Reads a partition from a binary file that contains an unsigned integer
	of the given width for each node.

	Parameters
	----------
	width : int, optional
		The width of the unsigned integer in bytes (4 or 8). Default: 4
	"""
	cdef _BinaryPartitionReader _this

	def __cinit__(self, uint8_t width=4):
		self._this = _BinaryPartitionReader(width)

	def read(self, path):
		"""
		read(path)
		
		Reads a partition from a binary file.
		
		Parameters
		----------
		path : str
			The input path.

		Returns
		-------
		networkit.Partition
			The resulting partition.
		"""
		return Partition().setThis(self._this.read(stdstring(path)))

cdef extern from "<networkit/io/BinaryPartitionWriter.hpp>":

	cdef cppclass _BinaryPartitionWriter "NetworKit::BinaryPartitionWriter":
		_BinaryPartitionWriter() except +
		_BinaryPartitionWriter(uint8_t width) except +
		_Partition write(_Partition zeta, string path) nogil except +

cdef class BinaryPartitionWriter:
	"""
	BinaryPartitionWriter(width=4)

	Writes a partition to a file to contains a binary list of partition ids.
	Partition ids are unsigned integers.

	Parameters
	----------
	width : int, optional
		The width of the unsigned integer in bytes (4 or 8). Default: 4

	"""
	cdef _BinaryPartitionWriter _this

	def __cinit__(self, uint8_t width=4):
		self._this = _BinaryPartitionWriter(width)

	def write(self, Partition P not None, path):
		"""
		Write the partition to the given file.

		Parameters:
		-----------
		path : str
			The output path
		"""
		cdef string c_path = stdstring(path)

		with nogil:
			self._this.write(P._this, c_path)

		return self

cdef extern from "<networkit/io/EdgeListPartitionReader.hpp>":

	cdef cppclass _EdgeListPartitionReader "NetworKit::EdgeListPartitionReader":
		_EdgeListPartitionReader() except +
		_EdgeListPartitionReader(node firstNode, char sepChar) except +
		_Partition read(string path) except +


cdef class EdgeListPartitionReader:
	""" 
	EdgeListPartitionReader(firstNode=1, sepChar = '`\`t')
	
	Reads a partition from an edge list type of file.

	Parameters
	----------
	firstNode : int, optional
		Id of first node. Default: 1
	sepChar : str
		Character which is used for data seperation. Default: '\t'
	"""
	cdef _EdgeListPartitionReader _this

	def __cinit__(self, node firstNode=1, sepChar = `'\t'`):
		self._this = _EdgeListPartitionReader(firstNode, stdstring(sepChar)[0])

	def read(self, path):
		"""
		read(path)
		
		Reads a partition from ad edge list file.
		
		Parameters
		----------
		path : str
			The input path.

		Returns
		-------
		networkit.Partition
			The resulting partition.
		"""
		return Partition().setThis(self._this.read(stdstring(path)))

cdef extern from "<networkit/io/BinaryEdgeListPartitionReader.hpp>":

	cdef cppclass _BinaryEdgeListPartitionReader "NetworKit::BinaryEdgeListPartitionReader":
		_BinaryEdgeListPartitionReader() except +
		_BinaryEdgeListPartitionReader(node firstNode, uint8_t width) except +
		_Partition read(string path) nogil except +
		_Partition read(vector[string] paths) nogil except +

cdef class BinaryEdgeListPartitionReader:
	"""
	BinaryEdgeListPartitionReader(firstNode=0, width=4)

	Reads a partition file that contains a binary list of pairs (node, partition(node)).
	It is assumed that all integers are unsigned.

	Parameters
	----------
	firstNode : int, optional
		The id of the first node, this is subtracted from all read node ids. Default: 0
	width : int, optional
		The width of the unsigned integer in bytes (4 or 8). Default: 4
	"""
	cdef _BinaryEdgeListPartitionReader _this

	def __cinit__(self, node firstNode=0, uint8_t width=4):
		self._this = _BinaryEdgeListPartitionReader(firstNode, width)

	def read(self, paths):
		"""
		read(paths)

		Read the partition from one or multiple files

		Parameters
		----------
		paths : str or list(str)
			The input path(s)

		Returns
		-------
		networkit.Partition
			The resulting partition.
		"""
		cdef vector[string] c_paths

		if isinstance(paths, str):
			c_paths.push_back(stdstring(paths))
		else:
			c_paths.reserve(len(paths))

			for p in paths:
				c_paths.push_back(stdstring(p))

		cdef _Partition result

		with nogil:
			result = move(self._this.read(c_paths)) # extra move in order to avoid copying the internal variable that is used by Cython

		return Partition().setThis(result)

cdef extern from "<networkit/io/BinaryEdgeListPartitionWriter.hpp>":

	cdef cppclass _BinaryEdgeListPartitionWriter "NetworKit::BinaryEdgeListPartitionWriter":
		_BinaryEdgeListPartitionWriter() except +
		_BinaryEdgeListPartitionWriter(node firstNode, uint8_t width) except +
		_Partition write(_Partition P, string path) nogil except +

cdef class BinaryEdgeListPartitionWriter:
	"""
	BinaryEdgeListPartitionWriter(firstNode=0, width=4)
	
	Writes a partition file that contains a binary list of pairs (node, partition(node)).

	Parameters
	----------
	firstNode : int, optional
		The id of the first node, this is added to all writen node ids. Default: 0
	width : int, optional
		The width of the unsigned integer in bytes (4 or 8). Default: 4
	"""
	cdef _BinaryEdgeListPartitionWriter _this

	def __cinit__(self, node firstNode=0, uint8_t width=4):
		self._this = _BinaryEdgeListPartitionWriter(firstNode, width)

	def write(self, Partition P not None, path):
		"""
		write(P, path)

		Write the partition to the given file.

		Parameters
		----------
		P : networkit.Partition
			The input partition.
		path : str
			The output path.
		"""
		cdef string c_path = stdstring(path)

		with nogil:
			self._this.write(P._this, c_path)

		return self

cdef extern from "<networkit/io/SNAPEdgeListPartitionReader.hpp>":

	cdef cppclass _SNAPEdgeListPartitionReader "NetworKit::SNAPEdgeListPartitionReader":
		_SNAPEdgeListPartitionReader() except +
		_Cover read(string path, unordered_map[node,node] nodeMap,_Graph G) except +
#		_Partition readWithInfo(string path, count nNodes) except +

cdef class SNAPEdgeListPartitionReader:
	""" 
	SNAPEdgeListPartitionReader(path, nodeMap, G)
	
	Reads a partition from a SNAP 'community with ground truth' file
	"""
	cdef _SNAPEdgeListPartitionReader _this

	def read(self,path, nodeMap, Graph G):
		cdef unordered_map[node,node] cNodeMap
		for (key,val) in nodeMap:
			cNodeMap[key] = val
		return Cover().setThis(self._this.read(stdstring(path), cNodeMap, G._this))

cdef extern from "<networkit/io/CoverReader.hpp>":

	cdef cppclass _CoverReader "NetworKit::CoverReader":
		_CoverReader() except +
		_Cover read(string path,_Graph G) except +

cdef class CoverReader:
	""" 
	CoverReader(path, G)

	Reads a cover from a file
	File format: each line contains the space-separated node ids of a community

	Parameters
	----------
	path : str
		Input file path.
	G : networkit.Graph
		Graph corresponding to the community from path.

	Returns
	-------
	networkit.Cover
		The resulting cover of a graph.
	"""
	cdef _CoverReader _this

	def read(self, path, Graph G):
		return Cover().setThis(self._this.read(stdstring(path), G._this))

cdef extern from "<networkit/io/CoverWriter.hpp>":

	cdef cppclass _CoverWriter "NetworKit::CoverWriter":
		_CoverWriter() except +
		void write(_Cover, string path) nogil except +


cdef class CoverWriter:
	""" 
	CoverWriter(zeta, path)

	Writes a partition to a file.
	File format: each line contains the space-separated node ids of a community
	
	Parameters
	----------
	zeta : networkit.Partition
		The input partition.
	path : str
		The output path.
	"""
	cdef _CoverWriter _this

	def write(self, Cover zeta, path):
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(zeta._this, cpath)

cdef extern from "<networkit/io/EdgeListCoverReader.hpp>":

	cdef cppclass _EdgeListCoverReader "NetworKit::EdgeListCoverReader":
		_EdgeListCoverReader() except +
		_EdgeListCoverReader(node firstNode) except +
		_Cover read(string path, _Graph G) except +


cdef class EdgeListCoverReader:
	""" 
	EdgeListCoverReader(firstNode=1)

	Reads a cover from an edge list type of file.
	File format: each line starts with a node id and continues with a list of the communities the node belongs to.
	
	Parameters
	----------
	firstNode : int, optional
		Id of first node. Default: 1
	"""
	cdef _EdgeListCoverReader _this

	def __cinit__(self, firstNode=1):
		self._this = _EdgeListCoverReader(firstNode)

	def read(self, path, Graph G):
		"""
		read(path, G)
		
		Reads a cover from an edge list file.
		
		Parameters
		----------
		path : str
			The input path.
		G : networkit.Graph
			Graph corresponding to the community from path.

		Returns
		-------
		networkit.Cover
			Cover of graph.
		"""
		return Cover().setThis(self._this.read(stdstring(path), G._this))

class __AutoNumber(Enum):
	def __new__(cls):
		value = len(cls.__members__) + 1
		obj = object.__new__(cls)
		obj._value_ = value
		return obj


class Format(__AutoNumber):
	""" 
	Simple enumeration class to list supported file types. Possible values: 

	- networkit.graphio.Format.DOT
	- networkit.graphio.Format.EdgeList
	- networkit.graphio.Format.EdgeListCommaOne
	- networkit.graphio.Format.EdgeListSpaceZero
	- networkit.graphio.Format.EdgeListSpaceOne
	- networkit.graphio.Format.EdgeListTabZero
	- networkit.graphio.Format.EdgeListTabOne
	- networkit.graphio.Format.GraphML
	- networkit.graphio.Format.GraphToolBinary
	- networkit.graphio.Format.GraphViz
	- networkit.graphio.Format.GEXF
	- networkit.graphio.Format.GML
	- networkit.graphio.Format.KONEC
	- networkit.graphio.Format.LFR
	- networkit.graphio.Format.METIS
	- networkit.graphio.Format.NetworkitBinary
	- networkit.graphio.Format.SNAP
	
	"""
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
	NetworkitBinary = ()

# reading

def getReader(fileformat, *kargs, **kwargs):
	"""
	getReader(fileformat, *kargs, **kwargs)

	Returns reader based on input fileformat.

	Parameters
	----------
	fileformat : networkit.graphio.Format
		A supported file format.
	`*kargs` : tuple()
		Additional input parameter (depending on the file format).
	`**kwargs` : dict()
		Additional input parameter (depending on the file format).
	"""
	#define your [edgelist] reader here:
	readers = {
		Format.METIS:			METISGraphReader(),
		Format.GraphML:			GraphMLReader(),
		Format.GEXF:			GEXFReader(),
		Format.SNAP:			SNAPGraphReader(),
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
		Format.NetworkitBinary:         NetworkitBinaryReader()
	}

	# special case for custom Edge Lists
	if fileformat == Format.EdgeList:
		if "continuous" in kwargs and kwargs["continuous"] == False:
			kwargs["firstNode"] = 0
		reader = EdgeListReader(*kargs, **kwargs)

	else:
		if fileformat not in readers:
			raise Exception("unrecognized format/format not supported as input: {0}".format(fileformat))
		reader = readers[fileformat]#(**kwargs)

	return reader


def readGraph(path, fileformat, *kargs, **kwargs):
	"""
	readGraph(path, fileformat, *kargs, **kwargs)

 	Read graph file in various formats and return a graph.

	Parameters
	----------
	fileformat : networkit.graphio.Format
		A supported file format.
	`*kargs` : tuple()
		Additional input parameter (depending on the file format).
	`**kwargs` : dict()
		Additional input parameter (depending on the file format). In case of a custom edge list, 
		pass the generic Fromat.EdgeList accompanied by the defining paramaters as follows:
		:code:`separator, firstNode, commentPrefix, continuous, directed`. :code:`commentPrefix, 
		continuous=True and directed` are optional because of their default values. 
		:code:`firstNode` is not needed when :code:`continuous=True`.
	"""
	reader = getReader(fileformat, *kargs, **kwargs)

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

def readGraphs(dirPath, pattern, fileformat, some=None, exclude=None, **kwargs):
	"""
	readGraphs(dirPath, pattern, fileformat, some=None, exclude=None, **kwargs)

 	Read all graph files contained in a directory whose filename contains the pattern, return 
	a dictionary of name to Graph object.

	Parameters
	----------
	dirPath : str
		Path, which contains input graphs.
	pattern : str
		Unix-style string pattern for file selection.
	fileformat : networkit.graphio.Format
		A supported file format.
	some : int, optional
		Restrict number of graphs to be read. Default: None
	exclude : str, optional
		Unix-style string pattern for file exclusion. Default: None
	`**kwargs` : dict()
		Additional input parameter (depending on the file format). In case of a custom edge list, 
		pass the generic Fromat.EdgeList accompanied by the defining paramaters as follows:
		:code:`separator, firstNode, commentPrefix, continuous, directed`. :code:`commentPrefix, 
		continuous=True and directed` are optional because of their default values. 
		:code:`firstNode` is not needed when :code:`continuous=True`.
	"""
	graphs = {}
	graph_id = 0
	for root, dirs, files in os.walk(dirPath):
		for file in files:
			if fnmatch.fnmatch(file, pattern):
				if (exclude is None) or (not fnmatch.fnmatch(file, exclude)):
					G = readGraph(os.path.join(root, file), fileformat, **kwargs)
					graphs[graph_id] = G
					graph_id += 1
					if some:
						if len(graphs) == some:
							return graphs
	return graphs


class MatReader:
	"""
	MatReader(key='G')

	Matlab file reader.
	File format: Adjacency matrix in matlab file format.

	Parameters
	----------
	key : str, optional
		Key to identify graph. Default: 'G'
	"""
	def __init__(self, key = "G"):
		self.key = key

	def read(self, path):
		"""
		read(path)
		
		Reads a graph from a matlab file.
		
		Parameters
		----------
		path : str
			The input path.

		Returns
		-------
		networkit.Graph
			The resulting graph.
		"""
		return readMat(path, self.key)

def readMat(path, key="G"):
	""" 
	readMat(key='G')

	Reads a Graph from a matlab object file containing an adjacency matrix and returns a graph.

	Parameters
	----------
	key : str, optional
		Key to identify graph. Default: 'G'
	"""
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
	"""
	MatWriter(key='G')

	Matlab file writer.

	Parameters
	----------
	key : str, optional
		Key to identify graph. Default: 'G'
	"""
	def __init__(self, key="G"):
		self.key = key

	def write(self, G, path, key="G"):
		"""
		write(G, path, key='G')
		
		Writes a graph to a file.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		path : str
			The output path.
		key : str, optional
			Key to identify graph. Default: 'G'
		"""
		writeMat(G, path, key)

def writeMat(G, path, key="G"):
	"""
	writeMat(G, path, key='G')
	
	Writes a graph to a file.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	path : str
		The output path.
	key : str, optional
		Key to identify graph. Default: 'G'
	"""
	matrix = algebraic.adjacencyMatrix(G, matrixType='sparse')
	scipy.io.savemat(path, {key : matrix})


# writing
def getWriter(fileformat, *kargs, **kwargs):
	"""
	getWriter(fileformat, *kargs, **kwargs)

	Returns reader based on input fileformat.

	Parameters
	----------
	fileformat : networkit.graphio.Format
		A supported file format.
	`*kargs` : tuple()
		Additional input parameter (depending on the file format).
	`**kwargs` : dict()
		Additional input parameter (depending on the file format).
	"""	
	writers =	{
		Format.METIS:			METISGraphWriter(),
		Format.GraphML:			GraphMLWriter(),
		Format.GEXF:			GEXFWriter(),
		Format.SNAP:			SNAPGraphWriter(),
		Format.EdgeListCommaOne:	EdgeListWriter(',',1,),
		Format.EdgeListSpaceOne:	EdgeListWriter(' ',1),
		Format.EdgeListSpaceZero:	EdgeListWriter(' ',0),
		Format.EdgeListTabOne:		EdgeListWriter('\t',1),
		Format.EdgeListTabZero:		EdgeListWriter('\t',0),
		Format.GraphViz:		DotGraphWriter(),
		Format.DOT:			DotGraphWriter(),
		Format.GML:			GMLGraphWriter(),
		Format.LFR:			EdgeListWriter('\t',1),
		Format.GraphToolBinary:         GraphToolBinaryWriter(),
		Format.MAT:			MatWriter(),
		Format.ThrillBinary:		ThrillGraphBinaryWriter(),
		Format.NetworkitBinary:         NetworkitBinaryWriter()
	}

	# special case for custom Edge Lists
	if fileformat == Format.EdgeList:
		return EdgeListWriter(*kargs, **kwargs)

	else:
		if fileformat not in writers:
			raise Exception("format {0} currently not supported".format(fileformat))

		return writers[fileformat]#(**kwargs)

def writeGraph(G, path, fileformat, *kargs, **kwargs):
	""" 
	writeGraph(G, path, fileformat, *kargs, **kwargs)
	
	Write graph to various output formats.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	path : str
		Output file path.
	fileformat : networkit.graphio.Format
		A supported file format.
	`*kargs` : tuple()
		Additional input parameter (depending on the file format).
	`**kwargs` : dict()
		Additional input parameter (depending on the file format). In case of a custom edge list, 
		pass the generic Fromat.EdgeList accompanied by the defining paramaters as follows:
		:code:`separator, firstNode, commentPrefix, continuous, directed`. :code:`commentPrefix, 
		continuous=True and directed` are optional because of their default values. 
		:code:`firstNode` is not needed when :code:`continuous=True`.
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
	writer = getWriter(fileformat, *kargs, **kwargs)
	writer.write(G, path)
	logging.info("wrote graph {0} to file {1}".format(G, path))


class GraphConverter:
	"""
	GraphConverter(reader, writer)

	Converts a format input to another or the same file format. The execute the conversion
	call `convert`.

	Parameters
	----------
	reader : networkit.graphio.GraphReader
		A supported graph reader.
	writer : networkit.graphio.GraphWriter
		A supported graph writer.
	"""
	def __init__(self, reader, writer):
		self.reader = reader
		self.writer = writer

	def convert(self, inPath, outPath):
		"""
		convert(inPath, outPath)

		Execute the conversion.

		Parameters
		----------
		inPath : str
			The input path.
		outPath : str
			The output path.
		"""
		G = self.reader.read(inPath)
		self.writer.write(G, outPath)

	def __str__(self):
		return "GraphConverter: {0} => {0}".format(self.reader, self.writer)

def getConverter(fromFormat, toFormat):
	"""
	getConverter(fromFormat, toFormat)

	Returns a converter for a given set of file formats.

	Parameters
	----------
	fromFormat : networkit.graphio.Format
		Source for conversion.
	toFormat : networkit.graphio.Format
		Target for conversion.

	Returns
	-------
	networkit.graphio.GraphConverter
		Corresponding GraphConverter object.
	"""
	reader = getReader(fromFormat)
	writer = getWriter(toFormat)
	return GraphConverter(reader, writer)


def convertGraph(fromFormat, toFormat, fromPath, toPath=None):
	"""
	convertGraph(fromFormat, toFormat, fromPath, toPath=None)

	Converts a graph given by a set of file formats and path.

	Parameters
	----------
	fromFormat : networkit.graphio.Format
		Source for conversion.
	toFormat : networkit.graphio.Format
		Target for conversion.
	fromPath : str
		The input path.
	toPath : str, optional
		The output path. Default: None
	"""
	converter = getConverter(fromFormat, toFormat)
	if toPath is None:
		toPath = "{0}.{1}.graph".format(fromPath.split(".")[0], toFormat)
	converter.convert(fromPath, toPath)
	print("converted {0} to {1}".format(fromPath, toPath))

# dynamic

def readStream(path, mapped=True, baseIndex=0):
	"""
	readStream(path, mapped=True, baseIndex=0)	
	
	Read a graph event stream from a file.

	Parameters
	----------
	path : str
		The input path.
	mapped : bool, optional
		Indicates whether the ids should be mapped. Default: True
	baseIndex : int, optional
		Sets base index of nodes. Default: 0
	"""
	return DGSStreamParser(path, mapped, baseIndex).getStream()

def writeStream(stream, path):
	"""
	writeStream(stream, path)
	
	Write a graph event stream to a file.

	Parameters
	----------
	stream : networkit.dynamics.GraphEvent
		The input event stream.
	path : str
		The output path.
	"""
	DGSWriter().write(stream, path)
