# distutils: language=c++

from libc.stdint cimport uint8_t

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.map cimport map as cmap
from libcpp.unordered_map cimport unordered_map
from .helpers import stdstring

import os
import logging
import numpy
import scipy.io
import fnmatch
import queue
import xml.etree.cElementTree as ET
import xml.sax
from enum import Enum
from warnings import warn
from xml.dom import minidom

from .dynamics import DGSWriter, DGSStreamParser, GraphEvent
from .graph cimport _Graph, Graph
from .graph import Graph as __Graph
from .structures cimport _Cover, Cover, _Partition, Partition, count, index, node
from . import algebraic
from .support import MissingDependencyError

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Graph move( _Graph t ) nogil # specialized declaration as general declaration disables template argument deduction and doesn't work
	_Partition move( _Partition t) nogil

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
		cmap[string,node] getNodeMap() except +

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
		cdef cmap[string,node] cResult = (<_EdgeListReader*>(self._this)).getNodeMap()
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

class GEXFReader:
	"""
	GEXFReader()
	This class provides a function to read a file in the
	GEXF (Graph Exchange XML Format) format.
	For more details see: http://gexf.net/
	"""
	def __init__(self):
		""" Initializes the GEXFReader class """
		self.mapping = dict()
		self.g = Graph(0)
		self.weighted = False
		self.directed = False
		self.dynamic = False
		self.hasDynamicWeights = False
		self.q = queue.Queue()
		self.eventStream = []
		self.nInitialNodes = 0
		self.timeFormat = ""

	def read(self, fpath):
		""" 
		read(fpath)
		Reads and returns the graph object defined in fpath.
		
		Parameters
		----------
		fpath : str
			File path for GEXF-file.
		"""
		#0. Reset internal vars and parse the xml
		self.__init__()
		doc = minidom.parse(fpath)

		#1. Determine if graph is dynamic, directed and has dynamically changing weights
		graph = doc.getElementsByTagName("graph")[0]
		if (graph.getAttribute("defaultedgetype") == "directed"):
			self.directed = True
		if (graph.getAttribute("mode") == "dynamic"):
			self.dynamic = True
		if self.dynamic:
			self.timeFormat = graph.getAttribute("timeformat")
		attributes = graph.getElementsByTagName("attribute")
		for att in attributes:
			if att.getAttribute("id") == "weight":
				self.hasDynamicWeights = True
				self.weighted = True

		#2. Read nodes and map them to IDs defined in GEXF file
		nodes = doc.getElementsByTagName("node")
		for n in nodes:
			u = n.getAttribute("id")
			if self.dynamic:
				"""
				A GEXF ID can be a string. However, this version of parser accepts ids
				in only 2 formats:
				1. id = "0,1,2," etc.
				2. id = "n0, n1, n2" etc.
				So either an integer or an integer that has n prefix.
				Gephi generates its random graphs in 2nd format for example.
				"""
				_id = ""
				try:
					_id = int(u)
				except:
					_id = int(u[1:])
				# 2-way mapping to refer nodes back in mapDynamicNodes() method
				self.mapping[u] = _id
				self.mapping[_id] = u
				controlList = {'elementAdded': False, 'elementDeleted': False}
				spells = n.getElementsByTagName("spell")
				if len(spells) > 0:
					for s in spells:
						self.parseDynamics(s, "n", controlList, u)
				else:
					self.parseDynamics(n, "n", controlList, u)
			else:
				self.mapping[u] = self.nInitialNodes
				self.nInitialNodes +=1
		if self.dynamic:
			self.mapDynamicNodes()

		#3. Read edges and determine if graph is weighted
		edges = doc.getElementsByTagName("edge")
		for e in edges:
			u = e.getAttribute("source")
			v = e.getAttribute("target")
			w = "1.0"
			if e.hasAttribute("weight"):
				self.weighted = True
				w = e.getAttribute("weight")
			if self.dynamic:
				controlList = {'elementAdded': False, 'elementDeleted': False}
				spells = e.getElementsByTagName("spell")
				if len(spells) > 0:
					for s in spells:
						self.parseDynamics(s, "e", controlList, u, v, w)
				else:
					self.parseDynamics(e, "e", controlList, u, v, w)
			else:
				self.q.put((u, v, w))

		#4. Create graph object
		self.g = Graph(self.nInitialNodes, self.weighted, self.directed)

		#5. Add initial edges to the graph and sort the eventStream by time
		#5.1 Adding initial edges
		while not self.q.empty():
			edge = self.q.get()
			(u, v, w) = (edge[0], edge[1], float(edge[2]))
			self.g.addEdge(self.mapping[u], self.mapping[v], w)

		#5.2 Sorting the eventStream by time and adding timeStep between events that happen in different times
		self.eventStream.sort(key=lambda x:x[1])
		for i in range(1, len(self.eventStream)):
			if self.eventStream[i][1] != self.eventStream[i-1][1]:
				self.eventStream.append((GraphEvent(GraphEvent.TIME_STEP, 0, 0, 0), self.eventStream[i-1][1]))
		self.eventStream.sort(key=lambda x:x[1])
		self.eventStream = [event[0] for event in self.eventStream]
		return (self.g, self.eventStream)



	def parseDynamics(self, element, elementType, controlList,  u,  v = "0", w = "0"):
		"""
		parseDynamics(element, elementType, controlList,  u,  v = "0", w = "0")
		Determine the operations as follows:
		1. Element has start and not deleted before: Create add event
		2. Element has start and deleted before: Create restore event
		3. Element has end:Create del event
		4. If an element has end before start(or no start at all), add it to the initial graph
		5. For dynamic edges, simply go over the attvalues and create weight update events
		A dynamic element must be defined either using only spells
		or inline attributes. These 2 shouldn't be mixed.
		(For example, Gephi will treat them differently. It'll ignore the inline declaration
		if the same element also contains spells)
		Parameters
		----------
		element : str
			Element to add during reading a GEXF-file.
		elementType : str
			Element type ("n" for node or "e" for edge).
		controlList : dict
			Dict with elements indicate element properties. Example :code:`{'elementAdded': False, 'elementDeleted': False}`
		u : str
			Node u involved in element parsing.
		v : str, optional
			Node v involved in element parsing. Default: "0"
		w : str, optional
			Edgeweight w involved in element parsing. Default: "0"
		"""
		startTime = element.getAttribute("start")
		if startTime == "":
			startTime = element.getAttribute("startopen")
		endTime	= element.getAttribute("end")
		if endTime == "":
			endTime	= element.getAttribute("endopen")
		if self.timeFormat != "date":
			try:
				startTime = float(startTime)
			except:
				pass
			try:
				endTime = float(endTime)
			except:
				pass

		if startTime != "" and endTime != "":
			if startTime < endTime and not controlList['elementDeleted']:
				self.createEvent(startTime, "a"+elementType, u, v, w)
				controlList['elementAdded'] = True
			else:
				self.createEvent(startTime, "r"+elementType, u, v, w)
			self.createEvent(endTime, "d"+elementType, u, v, w)
			controlList['elementDeleted'] = True

		if startTime != "" and endTime == "":
			if controlList['elementDeleted']:
				self.createEvent(startTime, "r"+elementType, u, v, w)
			else:
				self.createEvent(startTime, "a"+elementType, u, v, w)
				controlList['elementAdded'] = True

	 	# Handle dynamic edge weights here
		if elementType == "e" and self.hasDynamicWeights:
			attvalues = element.getElementsByTagName("attvalue")
			# If a spell is traversed, attvalues are siblings
			if len(attvalues) == 0:
				attvalues = element.parentNode.parentNode.getElementsByTagName("attvalue")
			for att in attvalues:
				if att.getAttribute("for") == "weight":
					w = att.getAttribute("value")
					startTime = att.getAttribute("start")
					if startTime == "":
						startTime = att.getAttribute("startopen")
					if self.timeFormat != "date":
						startTime = float(startTime)
					# If this edge is not added, first weight update indicates edge addition
					if not controlList['elementAdded']:
						self.createEvent(startTime, "a"+elementType, u, v, w)
						controlList['elementAdded'] = True
					else:
						self.createEvent(startTime, "c"+elementType, u, v, w)

		if startTime == "":
			if not controlList['elementAdded']:
				if elementType == "n":
					self.mapping[u] = self.nInitialNodes
					self.nInitialNodes += 1
				else:
					self.q.put((u,v,w))
				controlList['elementAdded'] = True
			if endTime != "":
				self.createEvent(endTime, "d"+elementType, u, v, w)
				controlList['elementDeleted'] = True

	def createEvent(self, eventTime, eventType, u, v, w):
		"""
		createEvent(eventTime, eventType, u, v, w)
		Creates a NetworKit::GraphEvent from the supplied parameters
		and passes it to eventStream.
		Parameters
		----------
		eventTime : int
			Timestep indicating when the event happen (creating an order of events).
		eventType : str
			Abbreviation string representing a graph event. Should be one of the following:
			:code:`e, an, dn, rn, ae, re, de, ce`.
		u : int
			Id of node u involved in graph event.
		v : int
			Id of node v involved in graph event.
		w : float
			Edgeweight of edge between u and v.
		"""
		event, u = None, self.mapping[u]
		if eventType[1] == "e":
			v, w = self.mapping[v], float(w)
		if eventType == "an":
			event = GraphEvent(GraphEvent.NODE_ADDITION, u, 0, 0)
		elif eventType == "dn":
			event = GraphEvent(GraphEvent.NODE_REMOVAL, u, 0, 0)
		elif eventType == "rn":
			event = GraphEvent(GraphEvent.NODE_RESTORATION, u, 0, 0)
		elif eventType == "ae" or eventType == "re":
			event = GraphEvent(GraphEvent.EDGE_ADDITION, u, v, w)
		elif eventType == "de":
			event = GraphEvent(GraphEvent.EDGE_REMOVAL, u, v, w)
		elif eventType == "ce":
			event = GraphEvent(GraphEvent.EDGE_WEIGHT_UPDATE, u, v, w)
		self.eventStream.append((event, eventTime))

	def mapDynamicNodes(self):
		"""
		mapDynamicNodes()
		Node ID of a dynamic node must be determined before it's mapped to its GEXF ID.
		This requires processing the sorted eventStream and figuring out the addition order of the nodes.
		After that, node addition/deletion/restoration operations of this node must be readded to eventStream
		with correct mapping.
		Note
		----
		New mapping of a node can be equal to old mapping of a node. In order to prevent collisions,
		isMapped array must be maintained and controlled.
		"""
		nNodes = self.nInitialNodes
		nEvent = len(self.eventStream)
		isMapped = [False] * nEvent
		self.eventStream.sort(key=lambda x:x[1])
		for i in range(0, nEvent):
			event = self.eventStream[i]
			# Only the nodes with addition event will get remapped.
			if not isMapped[i] and event[0].type == GraphEvent.NODE_ADDITION:
				u = event[0].u
				self.mapping[self.mapping[u]] = nNodes
				# All the other events of that node comes after it's addition event
				for j in range(i, len(self.eventStream)):
					event = self.eventStream[j]
					if not isMapped[j] and event[0].u == u:
						mappedEvent = GraphEvent(event[0].type, self.mapping[self.mapping[u]], 0, 0)
						self.eventStream[j] = (mappedEvent, event[1])
						isMapped[j] = True
				nNodes +=1
				isMapped[i] = True

	def getNodeMap(self):
		""" 
		getNodeMap()
		
		Returns
		-------
		dict(int ``:`` int)
			Dictionary containing mapping from GEXF ID to node ID
		"""
		forwardMap = dict()
		for key in self.mapping:
			if type(key) == str:
				forwardMap[key] = self.mapping[key]
		return forwardMap


# GEXFWriter
class GEXFWriter:
	""" 
	GEXFWriter()
	
	This class provides a function to write a NetworKit graph to a file in the GEXF format. 
	"""

	def __init__(self):
		""" Initializes the class. """
		self.edgeIdctr = 0
		self.q = queue.Queue()
		self.hasDynamicWeight = False

	def write(self, graph, fname, eventStream = [], mapping = []):
		"""
		write(graph, fname, evenStream = [], mapping = [])
		Writes a graph to the specified file fname.
		
		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		fname : str 
			The desired file path and name to be written to.
		eventStream : list(networkit.dynamics.GraphEvent)
			Stream of events, each represented by networkit.dynamics.GraphEvent.
		mapping : list(int)
			Random node mapping.
		"""
		#0. Reset internal vars
		self.__init__()

		#1. Start with the root element and the right header information
		root = ET.Element('gexf')
		root.set("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance")
		root.set("xsi:schemaLocation","http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd")
		root.set('version', '1.2')

		#2. Create graph element with appropriate information
		graphElement = ET.SubElement(root,"graph")
		if graph.isDirected():
			graphElement.set('defaultedgetype', 'directed')
		else:
			graphElement.set('defaultedgetype', 'undirected')
		if len(eventStream) > 0:
			graphElement.set('mode', 'dynamic')
			graphElement.set('timeformat', 'double')
			for event in eventStream:
				if event.type == GraphEvent.EDGE_WEIGHT_UPDATE:
					dynamicAtt = ET.SubElement(graphElement, "attributes")
					dynamicAtt.set('class', 'edge')
					dynamicAtt.set('mode', 'dynamic')
					dynamicWeight = ET.SubElement(dynamicAtt, "attribute")
					dynamicWeight.set('id', 'weight')
					dynamicWeight.set('title', 'Weight')
					dynamicWeight.set('type', 'float')
					self.hasDynamicWeight = True
					break
		else:
			graphElement.set('mode', 'static')

		#3. Add nodes
		nodesElement = ET.SubElement(graphElement, "nodes")
		nNodes, idArray = 0, []
		#3.1 Count the # of nodes (inital + dynamic nodes)
		for event in eventStream:
			if event.type == GraphEvent.NODE_ADDITION:
				nNodes +=1
		nNodes += graph.numberOfNodes()
		for i in range(0, nNodes):
			idArray.append(i)
		# Optional:Map nodes to a random mapping if user provided one
		if (len(mapping) > 0):
			if(nNodes != len(mapping)):
				raise Exception('Size of nodes and mapping must match')
			else:
				for i in range(0, nNodes):
					idArray[i] = mapping[i]

		#3.2 Write nodes to the gexf file
		for n in range(nNodes):
			nodeElement = ET.SubElement(nodesElement,'node')
			nodeElement.set('id', str(idArray[n]))
			self.writeEvent(nodeElement, eventStream, n)

		#4. Add edges
		edgesElement = ET.SubElement(graphElement, "edges")
		#4.1 Put all edges into a queue(inital + dynamic edges)
		for u, v in graph.iterEdges():
			self.q.put((u, v, graph.weight(u, v)))
		for event in eventStream:
			if event.type == GraphEvent.EDGE_ADDITION:
				self.q.put((event.u, event.v, event.w))
		#4.2 Write edges to the gexf file
		while not self.q.empty():
			edgeElement = ET.SubElement(edgesElement,'edge')
			e = self.q.get()
			edgeElement.set('source', str(idArray[e[0]]))
			edgeElement.set('target', str(idArray[e[1]]))
			edgeElement.set('id', "{0}".format(self.edgeIdctr))
			self.edgeIdctr += 1
			if graph.isWeighted():
				edgeElement.set('weight', str(e[2]))
			self.writeEvent(edgeElement, eventStream, e)

		#5. Write the generated tree to the file
		tree = ET.ElementTree(root)
		tree.write(fname,"utf-8",True)

	def writeEvent(self, xmlElement, eventStream, graphElement):
		"""
		writeEvent(xmlElement, eventStream, graphElement)
		Write a single event. This is a supporting function and should normally not be called independently.
		Parameters
		----------
		xmlElement : xml.etree.cElementTree
			XML-encoded element, representing one GEXF-element.
		eventStream : list(networkit.dynamics.GraphEvent)
			Stream of events, each represented by networkit.dynamics.GraphEvent.
		graphElement : tuple(int, int, float)
			Tuple representing one graph element given by (node u, node v, edge weight w). 
		"""
		# A var that indicates if the event belongs the graph element we traverse on
		matched = False
		startEvents = [GraphEvent.NODE_ADDITION, GraphEvent.EDGE_ADDITION, GraphEvent.NODE_RESTORATION]
		endEvents = [GraphEvent.NODE_REMOVAL, GraphEvent.EDGE_REMOVAL]
		nodeEvents = [GraphEvent.NODE_ADDITION, GraphEvent.NODE_REMOVAL, GraphEvent.NODE_RESTORATION]
		edgeEvents = [GraphEvent.EDGE_ADDITION, GraphEvent.EDGE_REMOVAL, GraphEvent.EDGE_WEIGHT_UPDATE]
		spellTag, weightTag, operation = False, False, ""
		timeStep = 0
		spellsElement, attValuesElement = None, None

		for event in eventStream:
			if event.type == GraphEvent.TIME_STEP:
				timeStep += 1
			if type(graphElement) == type(0): #a node is an integer
				matched = (event.type in nodeEvents and event.u == graphElement)
			else:
				matched = (event.type in edgeEvents and (event.u == graphElement[0] and event.v == graphElement[1]))
			if matched:
				# Handle weight update seperately
				if event.type == GraphEvent.EDGE_WEIGHT_UPDATE:
					if not weightTag:
						attvaluesElement = ET.SubElement(xmlElement, "attvalues")
						weightTag = True
					attvalue = ET.SubElement(attvaluesElement, "attvalue")
					attvalue.set('for', 'weight')
					attvalue.set('value', str(event.w))
					attvalue.set('start', str(timeStep))
					attvalue.set('endopen', str(timeStep + 1))
				else:
					if event.type in startEvents:
						operation = "start"
					else:
						operation = "end"
					if not spellTag:
						spellsElement = ET.SubElement(xmlElement, "spells")
						spellTag = True
					spellElement = ET.SubElement(spellsElement, "spell")
					spellElement.set(operation, str(timeStep))

class GraphMLSAX(xml.sax.ContentHandler):
	""" 
	GraphMLSAX()

	Parser for GraphML XML files, based on Pythons XML.SAX implementation.
	"""

	def __init__(self):
		""" Initializes several important variables """
		xml.sax.ContentHandler.__init__(self)
		self.charBuffer = []
		self.mapping = dict()
		self.g = Graph(0)
		self.graphName = 'unnamed'
		self.weightedID = ''
		self.weighted = False
		self.directed = False
		self.edgestack = []
		self.edgeweight = 0.0
		self.keepData = False

	def startElement(self, name, attrs):
		""" 
		startElement(name, attrs)

		Parses all currently relevant XML tags and retrieves data.
		
		Parameters
		----------
		name : str
			Name of the element. Possible values: graph, node, edge, key, data
		attr : dict()
			Attributes of element.
		"""
		if name == "graph":
			# determine, if graph is directed:
			if attrs.getValue("edgedefault") == "directed":
				print("identified graph as directed")
				self.directed = True
			if "id" in  attrs.getNames() and not attrs.getValue("id") == '':
					self.graphName = attrs.getValue("id")
			self.g = Graph(0,self.weighted, self.directed)
		if name == "node":
			u = self.g.addNode()
			val = attrs.getValue("id")
			self.mapping[val] = u
		elif name == "edge":
			u = attrs.getValue("source")
			v = attrs.getValue("target")
			self.edgestack.append((u,v))
		elif name == "key":
			#print("found element with tag KEY")
			if (attrs.getValue("for") == 'edge' and attrs.getValue("attr.name") == 'weight' and attrs.getValue("attr.type") == 'double'):
				self.weighted = True
				self.weightedID = attrs.getValue("id")
				print("identified graph as weighted")
		elif name == "data" and attrs.getValue("key") == self.weightedID:
			self.keepData = True

	def endElement(self, name):
		""" 
		endElement(name)

		Finalizes parsing of the started Element and processes retrieved data.

		Parameters
		----------
		name : str
			Name of the element. Possible values: edge, data
		"""
		data = self.getCharacterData()
		if name == "edge":
			u = self.edgestack[len(self.edgestack)-1][0]
			v = self.edgestack[len(self.edgestack)-1][1]
			self.edgestack.pop()
			if self.weighted:
				#print ("identified edge as weighted with weight: {0}".format(edgeweight))
				self.g.addEdge(self.mapping[u], self.mapping[v], self.edgeweight)
				self.edgeweight = 0.0
			else:
				self.g.addEdge(self.mapping[u], self.mapping[v])
		elif name == "data" and self.keepData:
			self.keepData = False
			self.edgeweight = float(data)

	def characters(self, content):
		""" 
		characters(content)

		Appends content string to the textbuffer.

		Parameters
		----------
		content : str
			String to be added.
		"""
		self.charBuffer.append(content)

	def getCharacterData(self):
		""" 
		getCharacterData()

		Returns current textbuffer and clears it afterwards.
		"""
		data = ''.join(self.charBuffer).strip()
		self.charBuffer = []
		return data

	def getGraph(self):
		""" 
		getGraph()

		Return parsed graph.
		"""
		return self.g

class GraphMLReader:
	""" 
	GraphMLReader()

	This class serves as wrapper for the GraphMLSAX class
	which is able to parse a GraphML XML file and construct
	a graph.
	"""

	def __init__(self):
		""" Initializes the GraphMLSAX class """
		self.graphmlsax = GraphMLSAX()

	def read(self, fpath):
		""" 
		read(fpath)
		
		Parses a GraphML XML file and returns the constructed Graph
		
		Parameters
		----------
		fpath: str
			The path to the file as a String.
		"""
		xml.sax.parse(fpath, self.graphmlsax)
		return self.graphmlsax.getGraph()

# GraphMLWriter
class GraphMLWriter:
	""" 
	GraphMLWriter()

	This class provides a function to write a NetworKit graph to a file in the 
	GraphML format.
	"""
	
	def __init__(self):
		""" Initializes the class. """
		self.edgeIdCounter = 0
		self.dir_str = ''

	def write(self, graph, fname, nodeAttributes = {}, edgeAttributes = {}):
		""" 
		write(self, graph, fname, nodeAttributes = {}, edgeAttributes = {})
		
		Writes a NetworKit graph to the specified file fname. 
		
		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		fname: str
			The desired file path and name to be written to.
		nodeAttributes: dict(), optional
			Dictionary of node attributes in the form attribute name => list of attribute values. Default: {}
		edgeAttributes: dict(), optional
			Dictionary of edge attributes in the form attribute name => list of attribute values. Default: {}
		"""
		# reset some internal variables in case more graphs are written with the same instance
		self.edgeIdCounter = 0
		self.dir_str = ''

		if len(edgeAttributes) > 0 and not graph.hasEdgeIds():
			raise RuntimeError("Error, graph must have edge ids if edge attributes are given")

		# start with the root element and the right header information
		root = ET.Element('graphml')
		root.set("xmlnsi","http://graphml.graphdrawing.org/xmlns")
		root.set("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance")
		root.set("xsi:schemaLocation","http://graphml.graphdrawing.org/xmlns \
			http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd")

		maxAttrKey = 1
		# if the graph is weighted, add the attribute
		if graph.isWeighted():
			attrElement = ET.SubElement(root,'key')
			attrElement.set('for','edge')
			attrElement.set('id', 'd1')
			attrElement.set('attr.name','weight')
			attrElement.set('attr.type','double')
			maxAttrKey += 1

		attrKeys = {}
		import numbers
		import itertools
		for attType, attName, attData in itertools.chain(
			map(lambda d : ('node', d[0], d[1]), nodeAttributes.items()),
			map(lambda d : ('edge', d[0], d[1]), edgeAttributes.items())):

			attrElement = ET.SubElement(root, 'key')
			attrElement.set('for', attType)
			attrElement.set('id', 'd{0}'.format(maxAttrKey))
			attrKeys[(attType, attName)] = 'd{0}'.format(maxAttrKey)
			maxAttrKey += 1
			attrElement.set('attr.name', attName)
			if len(attData) == 0:
				attrElement.set('attr.type', 'int')
			elif isinstance(attData[0], bool):
				attrElement.set('attr.type', 'boolean')
				# special handling for boolean attributes: convert boolean into lowercase string
				if attType == 'edge':
					edgeAttributes[attName] = [ str(d).lower() for d in attData ]
				else:
					nodeAttributes[attName] = [ str(d).lower() for d in attData ]
			elif isinstance(attData[0], numbers.Integral):
				attrElement.set('attr.type', 'int')
			elif isinstance(attData[0], numbers.Real):
				attrElement.set('attr.type', 'double')
			else:
				attrElement.set('attr.type', 'string')


		# create graph element with appropriate information
		graphElement = ET.SubElement(root,"graph")
		if graph.isDirected():
			graphElement.set('edgedefault', 'directed')
			self.dir_str = 'true'
		else:
			graphElement.set('edgedefault', 'undirected')
			self.dir_str = 'false'

		# Add nodes
		for n in graph.iterNodes():
			nodeElement = ET.SubElement(graphElement,'node')
			nodeElement.set('id', str(n))
			for attName, attData in nodeAttributes.items():
				dataElement = ET.SubElement(nodeElement, 'data')
				dataElement.set('key', attrKeys[('node', attName)])
				dataElement.text = str(attData[n])

		# in the future: more attributes
	        #for a in n.attributes():
        	#    if a != 'label':
	        #        data = doc.createElement('data')
        	#        data.setAttribute('key', a)
	        #        data.appendChild(doc.createTextNode(str(n[a])))
        	#        node.appendChild(data)

		# Add edges
		def addEdge(u, v, w, eid):
			edgeElement = ET.SubElement(graphElement,'edge')
			edgeElement.set('directed', self.dir_str)
			edgeElement.set('target', str(v))
			edgeElement.set('source', str(u))
			if graph.hasEdgeIds():
				edgeElement.set('id', "e{0}".format(eid))
			else:
				edgeElement.set('id', "e{0}".format(self.edgeIdCounter))
				self.edgeIdCounter += 1
			if graph.isWeighted():
				# add edge weight
				dataElement = ET.SubElement(edgeElement,'data')
				dataElement.set('key','d1')
				dataElement.text = str(w)
			for attName, attData in edgeAttributes.items():
				dataElement = ET.SubElement(edgeElement, 'data')
				dataElement.set('key', attrKeys[('edge', attName)])
				dataElement.text = str(attData[eid])
		graph.forEdges(addEdge)

	#TODO: optional prettify function for formatted output of xml files
		tree = ET.ElementTree(root)
		tree.write(fname,"utf-8",True)
