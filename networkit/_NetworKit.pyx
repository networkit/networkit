# cython: language_level=3

#includes
# needed for collections.Iterable
from networkit.exceptions import ReducedFunctionalityWarning
import collections
import math
import os
import tempfile
import warnings

try:
	import pandas
except:
	warnings.warn("WARNING: module 'pandas' not found, some functionality will be restricted", ReducedFunctionalityWarning)

# C++ operators
from cython.operator import dereference, preincrement

# type imports
from libc.stdint cimport uint64_t
from libc.stdint cimport int64_t
from libc.stdint cimport uint8_t

# the C++ standard library
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.stack cimport stack
from libcpp.string cimport string
from libcpp.unordered_set cimport unordered_set
from libcpp.unordered_map cimport unordered_map
from libcpp.algorithm cimport sort as stdsort

# NetworKit typedefs
ctypedef uint64_t count
ctypedef uint64_t index
ctypedef uint64_t edgeid
ctypedef index node
ctypedef index cluster
ctypedef double edgeweight
ctypedef double coordinate

from .base cimport _Algorithm
from .base cimport Algorithm
from .graph cimport _Graph, Graph
from .structures cimport _Cover, Cover, _Partition, Partition
from .matching cimport _Matching, Matching
from networkit.graphtools import GraphTools
from .dynamics cimport _GraphEvent, GraphEvent
from .community cimport _CommunityDetectionAlgorithm, CommunityDetector

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

none = _none

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes

def pystring(stdstring):
	""" convert a std::string (= python byte string) to a normal Python string"""
	return stdstring.decode("utf-8")

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Graph move( _Graph t ) nogil # specialized declaration as general declaration disables template argument deduction and doesn't work
	_Partition move( _Partition t) nogil
	pair[_Graph, vector[node]] move(pair[_Graph, vector[node]]) nogil
	vector[pair[pair[node, node], double]] move(vector[pair[pair[node, node], double]]) nogil
	vector[double] move(vector[double])
	vector[bool_t] move(vector[bool_t])
	vector[pair[node, node]] move(vector[pair[node, node]]) nogil

cdef extern from "<networkit/auxiliary/Parallel.hpp>" namespace "Aux::Parallel":

	void sort[Iter](Iter begin, Iter end) nogil
	void sort[Iter, Comp](Iter begin, Iter end, Comp compare) nogil

# Function definitions

cdef extern from "<networkit/auxiliary/Log.hpp>" namespace "Aux":

	#void _configureLogging "Aux::configureLogging" (string loglevel)
	string _getLogLevel "Aux::Log::getLogLevel" () except +
	void _setLogLevel "Aux::Log::setLogLevel" (string loglevel) except +
	void _setPrintLocation "Aux::Log::Settings::setPrintLocation" (bool_t) except +

def getLogLevel():
	""" Get the current log level"""
	return pystring(_getLogLevel())

def setLogLevel(loglevel):
	""" Set the current loglevel"""
	_setLogLevel(stdstring(loglevel))

def setPrintLocation(flag):
	""" Switch locations in log statements on or off"""
	_setPrintLocation(flag)

cdef extern from "<networkit/auxiliary/Parallelism.hpp>" namespace "Aux":

	void _setNumberOfThreads "Aux::setNumberOfThreads" (int)
	int _getCurrentNumberOfThreads "Aux::getCurrentNumberOfThreads" ()
	int _getMaxNumberOfThreads "Aux::getMaxNumberOfThreads" ()
	void _enableNestedParallelism "Aux::enableNestedParallelism" ()

def setNumberOfThreads(nThreads):
	""" Set the number of OpenMP threads """
	_setNumberOfThreads(nThreads)

def getCurrentNumberOfThreads():
	""" Get the number of currently running threads"""
	return _getCurrentNumberOfThreads()

def getMaxNumberOfThreads():
	""" Get the maximum number of available threads"""
	return _getMaxNumberOfThreads()

def enableNestedParallelism():
	""" Enable nested parallelism for OpenMP"""
	from warnings import warn
	warn("Nested parallelism has been deprecated.")

cdef extern from "<networkit/auxiliary/Random.hpp>" namespace "Aux::Random":

	void _setSeed "Aux::Random::setSeed" (uint64_t, bool_t)

def setSeed(uint64_t seed, bool_t useThreadId):
	""" Set the random seed that is used in NetworKit.

	Note that there is a separate random number generator per thread.

	Parameters
	----------
	seed : uint64_t
		The seed
	useThreadId : bool
		If the thread id shall be added to the seed
	"""
	_setSeed(seed, useThreadId)

# Class definitions

## Module: engineering

# TODO: timer

cdef extern from "<networkit/viz/Point.hpp>" namespace "NetworKit" nogil:

	cdef cppclass Point[T]:
		Point()
		Point(T x, T y)
		T& operator[](const index i) except +
		T& at(const index i) except +

	cdef cppclass _Point2D "NetworKit::Point2D":
		_Point2D()
		pair[coordinate, coordinate] asPair()

cdef object toPoint2DVector(const vector[_Point2D]& v):
	return [v[i].asPair() for i in range(v.size())]

cdef object toNodePoint2DVector(const vector[pair[node, _Point2D]]& v):
	return [(v[i].first, v[i].second.asPair()) for i in range(v.size())]

cdef extern from "<networkit/independentset/Luby.hpp>":

	cdef cppclass _Luby "NetworKit::Luby":
		_Luby() except +
		vector[bool_t] run(_Graph G) except +
		string toString()


# FIXME: check correctness
cdef class Luby:
	""" Luby's parallel maximal independent set algorithm"""
	cdef _Luby _this

	def run(self, Graph G not None):
		""" Returns a bool vector of length n where vec[v] is True iff v is in the independent sets.
		Parameters
		----------
		G : networkit.Graph
			The graph.
		Returns
		-------
		vector
			A bool vector of length n.
		"""
		return self._this.run(G._this)
		# TODO: return self

	def toString(self):
		""" Get string representation of the algorithm.
		Returns
		-------
		string
			The string representation of the algorithm.
		"""
		return self._this.toString().decode("utf-8")

# Module: generators


# cdef extern from "<networkit/generators/MultiscaleGenerator.hpp>":

# 	cdef cppclass _MultiscaleGenerator "NetworKit::MultiscaleGenerator":
# 		_MultiscaleGenerator(_Graph O) except +
# 		_Graph generate() except +
#
#
# cdef class MultiscaleGenerator:
# 	"""
# 	TODO:
# 	"""
# 	cdef _MultiscaleGenerator *_this
# 	cdef Graph O	# store reference to input graph to not let it be garbage-collection
#
# 	def __cinit__(self, Graph O):
# 		self._this = new _MultiscaleGenerator(O._this)
# 		self.O = O
#
# 	def generate(self):
# 		return Graph(0).setThis(self._this.generate())
#
# 	@classmethod
# 	def fit(cls, Graph G):
# 		return cls(G)




# Module: graphio

cdef extern from "<networkit/io/GraphReader.hpp>":

	cdef cppclass _GraphReader "NetworKit::GraphReader":
		_GraphReader() nogil except +
		_Graph read(string path) nogil except +

cdef extern from "<networkit/io/GraphReader.hpp>" namespace "NetworKit::GraphReader":

	cdef enum _MultipleEdgesHandling "NetworKit::GraphReader::MultipleEdgesHandling":
		DISCARD_EDGES,
		SUM_WEIGHTS_UP,
		KEEP_MINIMUM_WEIGHT

class MultipleEdgesHandling:
	DiscardEdges = DISCARD_EDGES
	SumWeightsUp = SUM_WEIGHTS_UP
	KeepMinimumWeight = KEEP_MINIMUM_WEIGHT

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
		cdef string cpath = stdstring(path)
		cdef _Graph result

		with nogil:
			result = move(self._this.read(cpath)) # extra move in order to avoid copying the internal variable that is used by Cython
		return Graph(0).setThis(result)

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
			raise RuntimeError("Error, you may not use GraphReader directly, use a sub-class instead")

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
		G     : networkit.Graph
			The graph to write
		paths : str
			The output path
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
	""" Reads the METIS adjacency file format [1]. If the Fast reader fails,
		use readGraph(path, graphio.formats.metis) as an alternative.
		[1]: http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
	"""
	def __cinit__(self):
		self._this = new _METISGraphReader()

cdef extern from "<networkit/io/NetworkitBinaryReader.hpp>":
	cdef cppclass _NetworkitBinaryReader "NetworKit::NetworkitBinaryReader" (_GraphReader):
		_NetworkitBinaryReader() except +

cdef class NetworkitBinaryReader(GraphReader):
	"""
	Reads a graph written in the custom Networkit format documented in cpp/io/NetworkitGraph.md
	"""

	def __cinit__(self):
		self._this = new _NetworkitBinaryReader()

cdef extern from "<networkit/io/NetworkitBinaryWriter.hpp>":
	cdef cppclass _NetworkitBinaryWriter "NetworKit::NetworkitBinaryWriter" (_GraphWriter):
		_NetworkitBinaryWriter() except +

cdef class NetworkitBinaryWriter(GraphWriter):
	def __cinit__(self):
		self._this = new _NetworkitBinaryWriter()

cdef extern from "<networkit/io/GraphToolBinaryReader.hpp>":

	cdef cppclass _GraphToolBinaryReader "NetworKit::GraphToolBinaryReader" (_GraphReader):
		_GraphToolBinaryReader() except +

cdef class GraphToolBinaryReader(GraphReader):
	""" Reads the binary file format defined by graph-tool[1].
		[1]: http://graph-tool.skewed.de/static/doc/gt_format.html
	"""
	def __cinit__(self):
		self._this = new _GraphToolBinaryReader()

cdef extern from "<networkit/io/ThrillGraphBinaryReader.hpp>":

	cdef cppclass _ThrillGraphBinaryReader "NetworKit::ThrillGraphBinaryReader" (_GraphReader):
		_ThrillGraphBinaryReader(count n) except +
		_Graph read(vector[string] paths) nogil except +

cdef class ThrillGraphBinaryReader(GraphReader):
	"""
	Reads a graph format consisting of a serialized DIA of vector<uint32_t> from thrill.
	When the number of nodes is given, reading the graph is more efficient.
	Otherwise nodes are added to the graph as they are encountered.
	Edges must be present only in one direction.

	Parameters
	----------
	n : count
		The number of nodes
	"""
	def __cinit__(self, count n = 0):
		self._this = new _ThrillGraphBinaryReader(n)

	"""
	Read the graph from one or multiple files

	Parameters
	----------
	paths : str or list[str]
		The input path(s)
	"""
	def read(self, paths):
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
	""" Reads a graph from various text-based edge list formats.

	EdgeListReader(self, separator, firstNode, commentPrefix="#", continuous=True, directed=False)

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
	the smallest id used in the file.

	The file may also include line comments which start with the commentPrefix.

	Parameters
	----------
	separator : char
		The separator character. Must have length of exactly one.
	firstNode : node
		The id of the first node, this value will be subtracted from all node ids
	commentPrefix : string
		Lines starting with this prefix will be ignored
	continuous : bool
		File uses continuous node ids.
	directed : bool
		Treat input file as a directed graph.
	"""
	def __cinit__(self, separator, firstNode, commentPrefix="#", continuous=True, directed=False):
		if len(separator) != 1 or ord(separator[0]) > 255:
			raise RuntimeError("separator has to be exactly one ascii character");

		self._this = new _EdgeListReader(stdstring(separator)[0], firstNode, stdstring(commentPrefix), continuous, directed)

	def getNodeMap(self):
		""" Returns mapping of non-continuous files.

		The mapping is returned as dict (string -> node) projecting the original
		labels (as strings) to the reassigned integer node ids.
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
	""" Reader for the KONECT graph format, which is described in detail on the KONECT website[1].

		[1]: http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
	"""
	def __cinit__(self, remapNodes = False, handlingmethod = MultipleEdgesHandling.DiscardEdges):
		self._this = new _KONECTGraphReader(remapNodes, handlingmethod)

cdef extern from "<networkit/io/GMLGraphReader.hpp>":

	cdef cppclass _GMLGraphReader "NetworKit::GMLGraphReader"(_GraphReader):
		_GMLGraphReader() except +

cdef class GMLGraphReader(GraphReader):
	""" Reader for the GML graph format, which is documented here [1].

		[1]: http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf
 	"""
	def __cinit__(self):
		self._this = new _GMLGraphReader()

cdef extern from "<networkit/io/METISGraphWriter.hpp>":

	cdef cppclass _METISGraphWriter "NetworKit::METISGraphWriter" (_GraphWriter):
		_METISGraphWriter() except +


cdef class METISGraphWriter(GraphWriter):
	""" Writes graphs in the METIS format"""

	def __cinit__(self):
		self._this = new _METISGraphWriter()

cdef extern from "<networkit/io/GraphToolBinaryWriter.hpp>":

	cdef cppclass _GraphToolBinaryWriter "NetworKit::GraphToolBinaryWriter" (_GraphWriter):
		_GraphToolBinaryWriter() except +


cdef class GraphToolBinaryWriter(GraphWriter):
	""" Reads the binary file format defined by graph-tool[1].
		[1]: http://graph-tool.skewed.de/static/doc/gt_format.html
	"""
	def __cinit__(self):
		self._this = new _GraphToolBinaryWriter()

cdef extern from "<networkit/io/DotGraphWriter.hpp>":

	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter" (_GraphWriter):
		_DotGraphWriter() except +

cdef class DotGraphWriter(GraphWriter):
	""" Writes graphs in the .dot/GraphViz format"""
	def __cinit__(self):
		self._this = new _DotGraphWriter()

cdef extern from "<networkit/io/GMLGraphWriter.hpp>":

	cdef cppclass _GMLGraphWriter "NetworKit::GMLGraphWriter" (_GraphWriter):
		_GMLGraphWriter() except +


cdef class GMLGraphWriter(GraphWriter):
	""" Writes a graph and its coordinates as a GML file.[1]
		[1] http://svn.bigcat.unimaas.nl/pvplugins/GML/trunk/docs/gml-technical-report.pdf """

	def __cinit__(self):
		self._this = new _GMLGraphWriter()

cdef extern from "<networkit/io/EdgeListWriter.hpp>":

	cdef cppclass _EdgeListWriter "NetworKit::EdgeListWriter" (_GraphWriter):
		_EdgeListWriter() except +
		_EdgeListWriter(char separator, node firstNode, bool_t bothDirections) except +

cdef class EdgeListWriter(GraphWriter):
	""" Writes graphs in various edge list formats.

	Parameters
	----------
	separator : string
		The separator character.
	firstNode : node
		The id of the first node, this value will be added to all node ids
	bothDirections : bool, optional
		If undirected edges shall be written in both directions, i.e., as symmetric directed graph (default: false)
	"""

	def __cinit__(self, separator, firstNode, bool_t bothDirections = False):
		cdef char sep = stdstring(separator)[0]
		self._this = new _EdgeListWriter(sep, firstNode, bothDirections)


cdef extern from "<networkit/io/LineFileReader.hpp>":

	cdef cppclass _LineFileReader "NetworKit::LineFileReader":
		_LineFileReader() except +
		vector[string] read(string path)


cdef class LineFileReader:
	""" Reads a file and puts each line in a list of strings """
	cdef _LineFileReader _this

	def read(self, path):
		return self._this.read(stdstring(path))


cdef extern from "<networkit/io/SNAPGraphWriter.hpp>":
	cdef cppclass _SNAPGraphWriter "NetworKit::SNAPGraphWriter" (_GraphWriter):
		_SNAPGraphWriter() except +

cdef class SNAPGraphWriter(GraphWriter):
	""" Writes graphs in a format suitable for the Georgia Tech SNAP software [1]
		[1]: http://snap-graph.sourceforge.net/
	"""

	def __cinit__(self):
		self._this = new _SNAPGraphWriter()

cdef extern from "<networkit/io/SNAPGraphReader.hpp>":

	cdef cppclass _SNAPGraphReader "NetworKit::SNAPGraphReader"(_GraphReader):
		_SNAPGraphReader() except +
		_SNAPGraphReader(bool_t directed, bool_t remapNodes, count nodeCount)

cdef class SNAPGraphReader(GraphReader):
	""" Reads a graph from the SNAP graph data collection [1]
		[1]: http://snap.stanford.edu/data/index.html
	"""
	def __cinit__(self, directed = False, remapNodes = True, nodeCount = 0):
		self._this = new _SNAPGraphReader(directed, remapNodes, nodeCount)



cdef extern from "<networkit/io/PartitionReader.hpp>":

	cdef cppclass _PartitionReader "NetworKit::PartitionReader":
		_PartitionReader() except +
		_Partition read(string path) except +


cdef class PartitionReader:
	""" Reads a partition from a file.
		File format: line i contains subset id of element i.
	 """
	cdef _PartitionReader _this

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))


cdef extern from "<networkit/io/PartitionWriter.hpp>":

	cdef cppclass _PartitionWriter "NetworKit::PartitionWriter":
		_PartitionWriter() except +
		void write(_Partition, string path) nogil except +


cdef class PartitionWriter:
	""" Writes a partition to a file.
		File format: line i contains subset id of element i.
	 """
	cdef _PartitionWriter _this

	def write(self, Partition zeta, path):
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
	Reads a partition from a binary file that contains an unsigned integer
	of the given width for each node.

	Parameters
	----------
	width : int
		the width of the unsigned integer in bytes (4 or 8)

	"""
	cdef _BinaryPartitionReader _this

	def __cinit__(self, uint8_t width=4):
		self._this = _BinaryPartitionReader(width)

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))

cdef extern from "<networkit/io/BinaryPartitionWriter.hpp>":

	cdef cppclass _BinaryPartitionWriter "NetworKit::BinaryPartitionWriter":
		_BinaryPartitionWriter() except +
		_BinaryPartitionWriter(uint8_t width) except +
		_Partition write(_Partition zeta, string path) nogil except +

cdef class BinaryPartitionWriter:
	"""
	Writes a partition to a file to contains a binary list of partition ids.
	Partition ids are unsigned integers.

	Parameters
	----------
	width : int
		the width of the unsigned integer in bytes (4 or 8)

	"""
	cdef _BinaryPartitionWriter _this

	def __cinit__(self, uint8_t width=4):
		self._this = _BinaryPartitionWriter(width)

	def write(self, Partition P not None, path):
		"""
		Write the partition to the given file.

		Parameters
		----------
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
	""" Reads a partition from an edge list type of file
	 """
	cdef _EdgeListPartitionReader _this

	def __cinit__(self, node firstNode=1, sepChar = '\t'):
		self._this = _EdgeListPartitionReader(firstNode, stdstring(sepChar)[0])

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))

cdef extern from "<networkit/io/BinaryEdgeListPartitionReader.hpp>":

	cdef cppclass _BinaryEdgeListPartitionReader "NetworKit::BinaryEdgeListPartitionReader":
		_BinaryEdgeListPartitionReader() except +
		_BinaryEdgeListPartitionReader(node firstNode, uint8_t width) except +
		_Partition read(string path) nogil except +
		_Partition read(vector[string] paths) nogil except +

cdef class BinaryEdgeListPartitionReader:
	"""
	Reads a partition file that contains a binary list of pairs (node, partition(node)).
	It is assumed that all integers are unsigned.

	Parameters
	----------
	firstNode : node
		The id of the first node, this is subtracted from all read node ids
	width : int
		The width of the unsigned integer in bytes (4 or 8)
	"""
	cdef _BinaryEdgeListPartitionReader _this

	def __cinit__(self, node firstNode=0, uint8_t width=4):
		self._this = _BinaryEdgeListPartitionReader(firstNode, width)

	def read(self, paths):
		"""
		Read the partition from one or multiple files

		Parameters
		----------
		paths : str or list[str]
			The input path(s)
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
	Writes a partition file that contains a binary list of pairs (node, partition(node)).

	Parameters
	----------
	firstNode : node
		The id of the first node, this is added to all writen node ids
	width : int
		The width of the unsigned integer in bytes (4 or 8)
	"""
	cdef _BinaryEdgeListPartitionWriter _this

	def __cinit__(self, node firstNode=0, uint8_t width=4):
		self._this = _BinaryEdgeListPartitionWriter(firstNode, width)

	def write(self, Partition P not None, path):
		"""
		Write the partition to the given file.

		Parameters
		----------
		path : str
			The output path
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
	""" Reads a partition from a SNAP 'community with ground truth' file
	 """
	cdef _SNAPEdgeListPartitionReader _this

	def read(self,path, nodeMap, Graph G):
		cdef unordered_map[node,node] cNodeMap
		for (key,val) in nodeMap:
			cNodeMap[key] = val
		return Cover().setThis(self._this.read(stdstring(path), cNodeMap, G._this))

#	def readWithInfo(self,path,nNodes):
#		return Partition().setThis(self._this.readWithInfo(stdstring(path),nNodes))

#not existing yet, maybe in the future?
#cdef extern from "<networkit/io/EdgeListPartitionWriter.hpp>":

#	cdef cppclass _EdgeListPartitionWriter "NetworKit::EdgeListPartitionWriter":
#		_EdgeListPartitionWriter() except +
#		void write(_Partition, string path)


#cdef class EdgeListPartitionWriter:
#	""" Writes a partition to a edge list type of file.
#		File format: a line contains the element id and the subsed id of the element.
#	 """
#	cdef _EdgeListPartitionWriter _this

#	def Write(self, Partition zeta, path):
#		self._this.write(zeta._this, stdstring(path))

cdef extern from "<networkit/io/CoverReader.hpp>":

	cdef cppclass _CoverReader "NetworKit::CoverReader":
		_CoverReader() except +
		_Cover read(string path,_Graph G) except +

cdef class CoverReader:
	""" Reads a cover from a file
		File format: each line contains the space-separated node ids of a community
	 """
	cdef _CoverReader _this

	def read(self, path, Graph G):
		return Cover().setThis(self._this.read(stdstring(path), G._this))

cdef extern from "<networkit/io/CoverWriter.hpp>":

	cdef cppclass _CoverWriter "NetworKit::CoverWriter":
		_CoverWriter() except +
		void write(_Cover, string path) nogil except +


cdef class CoverWriter:
	""" Writes a partition to a file.
		File format: each line contains the space-separated node ids of a community
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
	""" Reads a cover from an edge list type of file
		File format: each line starts with a node id and continues with a list of the communities the node belongs to
	 """
	cdef _EdgeListCoverReader _this

	def __cinit__(self, firstNode=1):
		self._this = _EdgeListCoverReader(firstNode)

	def read(self, path, Graph G):
		return Cover().setThis(self._this.read(stdstring(path), G._this))



# Module: community

# Fused type for methods that accept both a partition and a cover
ctypedef fused PartitionCover:
	Partition
	Cover

cdef extern from "<networkit/community/ClusteringGenerator.hpp>":

	cdef cppclass _ClusteringGenerator "NetworKit::ClusteringGenerator":
		_ClusteringGenerator() except +
		_Partition makeSingletonClustering(_Graph G) except +
		_Partition makeOneClustering(_Graph G) except +
		_Partition makeRandomClustering(_Graph G, count k) except +
		_Partition makeContinuousBalancedClustering(_Graph G, count k) except +
		_Partition makeNoncontinuousBalancedClustering(_Graph G, count k) except +

cdef class ClusteringGenerator:
	""" Generators for various clusterings """
	cdef _ClusteringGenerator _this
	def makeSingletonClustering(self, Graph G):
		"""  Generate a clustering where each node has its own cluster

		Parameters
		----------
		G : networkit.Graph
			The graph for which the clustering shall be generated

		Returns
		-------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeSingletonClustering(G._this))
	def makeOneClustering(self, Graph G):
		"""  Generate a clustering with one cluster consisting of all nodes

		Parameters
		----------
		G : networkit.Graph
			The graph for which the clustering shall be generated

		Returns
		-------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeOneClustering(G._this))
	def makeRandomClustering(self, Graph G, count k):
		"""  Generate a clustering with `k` clusters to which nodes are assigned randomly

		Parameters
		----------
		G : networkit.Graph
			The graph for which the clustering shall be generated
		k: count
			The number of clusters that shall be generated

		Returns
		-------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeRandomClustering(G._this, k))
	def makeContinuousBalancedClustering(self, Graph G, count k):
		"""  Generate a clustering with `k` clusters to which nodes are assigned in continuous blocks

		Parameters
		----------
		G : networkit.Graph
			The graph for which the clustering shall be generated
		k: count
			The number of clusters that shall be generated

		Returns
		-------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeContinuousBalancedClustering(G._this, k))
	def makeNoncontinuousBalancedClustering(self, Graph G, count k):
		"""  Generate a clustering with `k` clusters, the ith node is assigned to cluster i % k. This means that
		for k**2 nodes, this clustering is complementary to the continuous clustering in the sense that no pair
		of nodes that is in the same cluster in one of the clusterings is in the same cluster in the other clustering.

		Parameters
		----------
		G : networkit.Graph
			The graph for which the clustering shall be generated
		k: count
			The number of clusters that shall be generated

		Returns
		-------
		networkit.Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeNoncontinuousBalancedClustering(G._this, k))

cdef extern from "<networkit/community/GraphClusteringTools.hpp>" namespace "NetworKit::GraphClusteringTools":

	float getImbalance(_Partition zeta) except +
	_Graph communicationGraph(_Graph graph, _Partition zeta) except +
	count weightedDegreeWithCluster(_Graph graph, _Partition zeta, node u, index cid)
	bool_t isProperClustering(_Graph G, _Partition zeta)
	bool_t isSingletonClustering(_Graph G, _Partition zeta)
	bool_t isOneClustering(_Graph G, _Partition zeta)
	bool_t equalClusterings(_Partition zeta, _Partition eta, _Graph G)

cdef class GraphClusteringTools:
	@staticmethod
	def getImbalance(Partition zeta):
		return getImbalance(zeta._this)
	@staticmethod
	def communicationGraph(Graph graph, Partition zeta):
		return Graph().setThis(communicationGraph(graph._this, zeta._this))
	@staticmethod
	def weightedDegreeWithCluster(Graph graph, Partition zeta, node u, index cid):
		return weightedDegreeWithCluster(graph._this, zeta._this, u, cid)
	@staticmethod
	def isProperClustering(Graph G, Partition zeta):
		return isProperClustering(G._this, zeta._this)
	@staticmethod
	def isSingletonClustering(Graph G, Partition zeta):
		return isSingletonClustering(G._this, zeta._this)
	@staticmethod
	def isOneClustering(Graph G, Partition zeta):
		return isOneClustering(G._this, zeta._this)
	@staticmethod
	def equalClustering(Partition zeta, Partition eta, Graph G):
		return equalClusterings(zeta._this, eta._this, G._this)


cdef extern from "<networkit/community/PartitionIntersection.hpp>":

	cdef cppclass _PartitionIntersection "NetworKit::PartitionIntersection":
		_PartitionIntersection() except +
		_Partition calculate(_Partition zeta, _Partition eta) except +

cdef class PartitionIntersection:
	""" Class for calculating the intersection of two partitions, i.e. the clustering with the fewest clusters
	such that each cluster is a subset of a cluster in both partitions.
	"""
	cdef _PartitionIntersection _this
	def calculate(self, Partition zeta, Partition eta):
		"""  Calculate the intersection of two partitions `zeta` and `eta`

		Parameters
		----------
		zeta: networkit.Partition
			The first partition
		eta: networkit.Partition
			The second partition

		Returns
		-------
		networkit.Partition
			The intersection of zeta and eta
		"""
		return Partition().setThis(self._this.calculate(zeta._this, eta._this))

cdef extern from "<networkit/community/Coverage.hpp>":

	cdef cppclass _Coverage "NetworKit::Coverage":
		_Coverage() except +
		double getQuality(_Partition _zeta, _Graph _G) except +

cdef class Coverage:
	""" Coverage is the fraction of intra-community edges """
	cdef _Coverage _this

	def getQuality(self, Partition zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)


cdef extern from "<networkit/community/EdgeCut.hpp>":

	cdef cppclass _EdgeCut "NetworKit::EdgeCut":
		_EdgeCut() except +
		double getQuality(_Partition _zeta, _Graph _G) except +

cdef class EdgeCut:
	""" Edge cut is the total weight of inter-community edges"""
	cdef _EdgeCut _this

	def getQuality(self, Partition zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)


cdef extern from "<networkit/community/Modularity.hpp>":

	cdef cppclass _Modularity "NetworKit::Modularity":
		_Modularity() except +
		double getQuality(_Partition _zeta, _Graph _G) nogil except +


cdef class Modularity:
	"""	Modularity is a quality index for community detection.
	It assigns a quality value in [-0.5, 1.0] to a partition of a graph which is higher for more modular networks and
	partitions which better capture the modular structure. See also http://en.wikipedia.org/wiki/Modularity_(networks).

 	Notes
	-----
	Modularity is defined as:

	.. math:: mod(\zeta) := \\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)} - \\frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }

	"""
	cdef _Modularity _this

	def getQuality(self, Partition zeta, Graph G):
		cdef double ret
		with nogil:
			ret = self._this.getQuality(zeta._this, G._this)
		return ret

cdef extern from "<networkit/community/HubDominance.hpp>":

	cdef cppclass _HubDominance "NetworKit::HubDominance":
		_HubDominance() except +
		double getQuality(_Partition _zeta, _Graph _G) except +
		double getQuality(_Cover _zeta, _Graph _G) except +

cdef class HubDominance:
	"""
	A quality measure that measures the dominance of hubs in clusters. The hub dominance of a single
	cluster is defined as the maximum cluster-internal degree of a node in that cluster divided by
	the maximum cluster-internal degree, i.e. the number of nodes in the cluster minus one. The
	value for all clusters is defined as the average of all clusters.

	Strictly speaking this is not a quality measure as this is rather dependent on the type of the
	considered graph, for more information see
	Lancichinetti A, Kivel M, Saramki J, Fortunato S (2010)
	Characterizing the Community Structure of Complex Networks
	PLoS ONE 5(8): e11976. doi: 10.1371/journal.pone.0011976
	http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011976
	"""

	cdef _HubDominance _this

	def getQuality(self, PartitionCover zeta, Graph G):
		"""
		Calculates the dominance of hubs in the given Partition or Cover of the given
		Graph.

		Parameters
		----------
		zeta : networkit.Partition or networkit.Cover
			The Partition or Cover for which the hub dominance shall be calculated
		G : networkit.Graph
			The Graph to which zeta belongs

		Returns
		-------
		double
			The average hub dominance in the given Partition or Cover
		"""
		return self._this.getQuality(zeta._this, G._this)

cdef extern from "<networkit/community/PLP.hpp>":

	cdef cppclass _PLP "NetworKit::PLP"(_CommunityDetectionAlgorithm):
		_PLP(_Graph _G, count updateThreshold, count maxIterations) except +
		_PLP(_Graph _G, _Partition baseClustering, count updateThreshold) except +
		count numberOfIterations() except +
		vector[count] getTiming() except +


cdef class PLP(CommunityDetector):
	""" Parallel label propagation for community detection:
	Moderate solution quality, very short time to solution.

	Parameters
	----------
	G : networkit.Graph
		The graph on which the algorithm has to run.
	updateThreshold : int
		number of nodes that have to be changed in each iteration so that a new iteration starts.
	baseClustering : networkit.Partition
		PLP needs a base clustering to start from; if none is given the algorithm will run on a singleton clustering.

	Notes
	-----
	As described in Ovelgoenne et al: An Ensemble Learning Strategy for Graph Clustering
 	Raghavan et al. proposed a label propagation algorithm for graph clustering.
 	This algorithm initializes every vertex of a graph with a unique label. Then, in iterative
 	sweeps over the set of vertices the vertex labels are updated. A vertex gets the label
 	that the maximum number of its neighbors have. The procedure is stopped when every vertex
 	has the label that at least half of its neighbors have.
	"""

	def __cinit__(self, Graph G not None, count updateThreshold=none, count maxIterations=none, Partition baseClustering=None,):
		"""
		Constructor to the Parallel label propagation community detection algorithm.

		"""
		self._G = G


		if baseClustering is None:
			self._this = new _PLP(G._this, updateThreshold, maxIterations)
		else:
			self._this = new _PLP(G._this, baseClustering._this, updateThreshold)


	def numberOfIterations(self):
		""" Get number of iterations in last run.

		Returns
		-------
		count
			The number of iterations.
		"""
		return (<_PLP*>(self._this)).numberOfIterations()

	def getTiming(self):
		""" Get list of running times for each iteration.

		Returns
		-------
		count
			The list of running times in milliseconds.
		"""
		return (<_PLP*>(self._this)).getTiming()

cdef extern from "<networkit/community/LPDegreeOrdered.hpp>":

	cdef cppclass _LPDegreeOrdered "NetworKit::LPDegreeOrdered"(_CommunityDetectionAlgorithm):
		_LPDegreeOrdered(_Graph _G) except +
		count numberOfIterations()

cdef class LPDegreeOrdered(CommunityDetector):
	""" Label propagation-based community detection algorithm which processes nodes in increasing order of node degree.	"""

	def __cinit__(self, Graph G not None):
		self._G = G
		self._this = new _LPDegreeOrdered(G._this)

	def numberOfIterations(self):
		""" Get number of iterations in last run.

		Returns
		-------
		count
			Number of iterations.
		"""
		return (<_LPDegreeOrdered*>(self._this)).numberOfIterations()

cdef extern from "<networkit/community/CutClustering.hpp>":

	cdef cppclass _CutClustering "NetworKit::CutClustering"(_CommunityDetectionAlgorithm):
		_CutClustering(_Graph _G) except +
		_CutClustering(_Graph _G, edgeweight alpha) except +

cdef extern from "<networkit/community/CutClustering.hpp>" namespace "NetworKit::CutClustering":

	map[double, _Partition] CutClustering_getClusterHierarchy "NetworKit::CutClustering::getClusterHierarchy"(const _Graph& G) nogil except +


cdef class CutClustering(CommunityDetector):
	"""
	Cut clustering algorithm as defined in
	Flake, Gary William; Tarjan, Robert E.; Tsioutsiouliklis, Kostas. Graph Clustering and Minimum Cut Trees.
	Internet Mathematics 1 (2003), no. 4, 385--408.

	Parameters
	----------
	G : networkit.Graph
	alpha : double
		The parameter for the cut clustering algorithm
	"""
	def __cinit__(self, Graph G not None,  edgeweight alpha):
		self._G = G
		self._this = new _CutClustering(G._this, alpha)

	@staticmethod
	def getClusterHierarchy(Graph G not None):
		""" Get the complete hierarchy with all possible parameter values.

		Each reported parameter value is the lower bound for the range in which the corresponding clustering is calculated by the cut clustering algorithm.

		Warning: all reported parameter values are slightly too high in order to avoid wrong clusterings because of numerical inaccuracies.
		Furthermore the completeness of the hierarchy cannot be guaranteed because of these inaccuracies.
		This implementation hasn't been optimized for performance.

		Parameters
		----------
		G : networkit.Graph
			The graph.

		Returns
		-------
		dict
			A dictionary with the parameter values as keys and the corresponding Partition instances as values
		"""
		cdef map[double, _Partition] result
		# FIXME: this probably copies the whole hierarchy because of exception handling, using move might fix this
		with nogil:
			result = CutClustering_getClusterHierarchy(G._this)
		pyResult = {}
		# FIXME: this code copies the partitions a lot!
		for res in result:
			pyResult[res.first] = Partition().setThis(res.second)
		return pyResult

cdef class DissimilarityMeasure:
	""" Abstract base class for partition/community dissimilarity measures """
	# TODO: use conventional class design of parametrized constructor, run-method and getters
	pass


cdef extern from "<networkit/community/NodeStructuralRandMeasure.hpp>":

	cdef cppclass _NodeStructuralRandMeasure "NetworKit::NodeStructuralRandMeasure":
		_NodeStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class NodeStructuralRandMeasure(DissimilarityMeasure):
	""" The node-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering all pairs of nodes.
	"""
	cdef _NodeStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret


cdef extern from "<networkit/community/GraphStructuralRandMeasure.hpp>":

	cdef cppclass _GraphStructuralRandMeasure "NetworKit::GraphStructuralRandMeasure":
		_GraphStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class GraphStructuralRandMeasure(DissimilarityMeasure):
	""" The graph-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering connected pairs of nodes.
	"""
	cdef _GraphStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret


cdef extern from "<networkit/community/JaccardMeasure.hpp>":

	cdef cppclass _JaccardMeasure "NetworKit::JaccardMeasure":
		_JaccardMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class JaccardMeasure(DissimilarityMeasure):
	""" TODO:
	"""
	cdef _JaccardMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret

cdef extern from "<networkit/community/NMIDistance.hpp>":

	cdef cppclass _NMIDistance "NetworKit::NMIDistance":
		_NMIDistance() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class NMIDistance(DissimilarityMeasure):
	""" The NMI distance assigns a similarity value in [0,1] to two partitions
		of a graph.
	"""
	cdef _NMIDistance _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret

cdef extern from "<networkit/community/AdjustedRandMeasure.hpp>":

	cdef cppclass _AdjustedRandMeasure "NetworKit::AdjustedRandMeasure":
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class AdjustedRandMeasure(DissimilarityMeasure):
	"""
	The adjusted rand dissimilarity measure as proposed by Huber and Arabie in "Comparing partitions" (http://link.springer.com/article/10.1007/BF01908075)
	"""
	cdef _AdjustedRandMeasure _this

	def getDissimilarity(self, Graph G not None, Partition first not None, Partition second not None):
		"""
		Get the adjust rand dissimilarity. Runs in O(n log(n)).

		Note that the dissimilarity can be larger than 1 if the partitions are more different than expected in the random model.

		Parameters
		----------
		G : networkit.Graph
			The graph on which the partitions shall be compared
		zeta : networkit.Partition
			The first partiton
		eta : networkit.Partition
			The second partition

		Returns
		-------
		double
			The adjusted rand dissimilarity
		"""
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret

cdef extern from "<networkit/community/LocalCommunityEvaluation.hpp>":

	cdef cppclass _LocalCommunityEvaluation "NetworKit::LocalCommunityEvaluation"(_Algorithm):
		double getWeightedAverage() except +
		double getUnweightedAverage() except +
		double getMaximumValue() except +
		double getMinimumValue() except +
		double getValue(index i) except +
		vector[double] getValues() except +
		bool_t isSmallBetter() except +

cdef class LocalCommunityEvaluation(Algorithm):
	"""
	Virtual base class of all evaluation methods for a single clustering which is based on the evaluation of single clusters.
	This is the base class both for Partitions as well as for Covers.
	"""
	def __init__(self, *args, **namedargs):
		if type(self) == LocalCommunityEvaluation:
			raise RuntimeError("Error, you may not use LocalCommunityEvaluation directly, use a sub-class instead")

	def getWeightedAverage(self):
		""" Get the average value weighted by cluster size.

		Returns
		-------
		double:
			The weighted average value.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getWeightedAverage()

	def getUnweightedAverage(self):
		""" Get the (unweighted) average value of all clusters.

		Returns
		-------
		double:
			The unweighted average value.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getUnweightedAverage()

	def getMaximumValue(self):
		""" Get the maximum value of all clusters.

		Returns
		-------
		double:
			The maximum value.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getMaximumValue()

	def getMinimumValue(self):
		""" Get the minimum value of all clusters.

		Returns
		-------
		double:
			The minimum value.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getMinimumValue()

	def getValue(self, index i):
		""" Get the value of the specified cluster.

		Parameters
		----------
		i : index
			The cluster to get the value for.

		Returns
		-------
		double:
			The value of cluster i.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getValue(i)

	def getValues(self):
		""" Get the values of all clusters.

		Returns
		-------
		list[double]:
			The values of all clusters.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).getValues()

	def isSmallBetter(self):
		""" If small values are better (otherwise large values are better).

		Returns
		-------
		bool:
			If small values are better.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_LocalCommunityEvaluation*>(self._this)).isSmallBetter()

cdef extern from "<networkit/community/LocalPartitionEvaluation.hpp>":

	cdef cppclass _LocalPartitionEvaluation "NetworKit::LocalPartitionEvaluation"(_LocalCommunityEvaluation):
		pass

cdef class LocalPartitionEvaluation(LocalCommunityEvaluation):
	"""
	Virtual base class of all evaluation methods for a single clustering which is based on the evaluation of single clusters.
	This is the base class for Partitions.
	"""
	cdef Graph _G
	cdef Partition _P

	def __init__(self, *args, **namedargs):
		if type(self) == LocalPartitionEvaluation:
			raise RuntimeError("Error, you may not use LocalPartitionEvaluation directly, use a sub-class instead")

	def __cinit__(self, Graph G not None, Partition P not None, *args, **namedargs):
		self._G = G
		self._P = P

	def __dealloc__(self):
		# Just to be sure that everything is properly deleted
		self._G = None
		self._P = None


cdef extern from "<networkit/community/LocalCoverEvaluation.hpp>":

	cdef cppclass _LocalCoverEvaluation "NetworKit::LocalCoverEvaluation"(_LocalCommunityEvaluation):
		pass


cdef class LocalCoverEvaluation(LocalCommunityEvaluation):
	"""
	Virtual base class of all evaluation methods for a single clustering which is based on the evaluation of single clusters.
	This is the base class for Covers.
	"""
	cdef Graph _G
	cdef Cover _C

	def __init__(self, *args, **namedargs):
		if type(self) == LocalCoverEvaluation:
			raise RuntimeError("Error, you may not use LocalCoverEvaluation directly, use a sub-class instead")

	def __cinit__(self, Graph G not None, Cover C not None, *args, **namedargs):
		self._G = G
		self._C = C

	def __dealloc__(self):
		# Just to be sure that everything is properly deleted
		self._G = None
		self._C = None

cdef extern from "<networkit/community/IntrapartitionDensity.hpp>":

	cdef cppclass _IntrapartitionDensity "NetworKit::IntrapartitionDensity"(_LocalPartitionEvaluation):
		_IntrapartitionDensity(_Graph G, _Partition P) except +
		double getGlobal() except +

cdef class IntrapartitionDensity(LocalPartitionEvaluation):
	"""
	The intra-cluster density of a partition is defined as the number of existing edges divided by the number of possible edges.
	The global value is the sum of all existing intra-cluster edges divided by the sum of all possible intra-cluster edges.

	Parameters
	----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _IntrapartitionDensity(self._G._this, self._P._this)

	def getGlobal(self):
		""" Get the global intra-cluster density.

		Returns
		-------
		double:
			The global intra-cluster density.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_IntrapartitionDensity*>(self._this)).getGlobal()


cdef extern from "<networkit/community/IsolatedInterpartitionConductance.hpp>":

	cdef cppclass _IsolatedInterpartitionConductance "NetworKit::IsolatedInterpartitionConductance"(_LocalPartitionEvaluation):
		_IsolatedInterpartitionConductance(_Graph G, _Partition P) except +

cdef class IsolatedInterpartitionConductance(LocalPartitionEvaluation):
	"""
	Isolated inter-partition conductance is a measure for how well a partition
	(communtiy/cluster) is separated from the rest of the graph.

	The conductance of a partition is defined as the weight of the cut divided
	by the volume (the sum of the degrees) of the nodes in the partition or the
	nodes in the rest of the graph, whatever is smaller. Small values thus indicate
	that the cut is small compared to the volume of the smaller of the separated
	parts. For the whole partitions usually the maximum or the unweighted average
	is used.

	See also Experiments on Density-Constrained Graph Clustering,
	Robert Grke, Andrea Kappes and  Dorothea Wagner, JEA 2015:
	http://dx.doi.org/10.1145/2638551

	Parameters
	----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _IsolatedInterpartitionConductance(self._G._this, self._P._this)

cdef extern from "<networkit/community/IsolatedInterpartitionExpansion.hpp>":

	cdef cppclass _IsolatedInterpartitionExpansion "NetworKit::IsolatedInterpartitionExpansion"(_LocalPartitionEvaluation):
		_IsolatedInterpartitionExpansion(_Graph G, _Partition P) except +

cdef class IsolatedInterpartitionExpansion(LocalPartitionEvaluation):
	"""
	Isolated inter-partition expansion is a measure for how well a partition
	(communtiy/cluster) is separated from the rest of the graph.

	The expansion of a partition is defined as the weight of the cut divided
	by number of nodes in the partition or in the rest of the graph, whatever
	is smaller. Small values thus indicate that the cut is small compared to
	the size of the smaller of the separated parts. For the whole partitions
	usually the maximum or the unweighted average is used. Note that expansion
	values can be larger than 1.

	See also Experiments on Density-Constrained Graph Clustering,
	Robert Grke, Andrea Kappes and Dorothea Wagner, JEA 2015:
	http://dx.doi.org/10.1145/2638551

	Parameters
	----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _IsolatedInterpartitionExpansion(self._G._this, self._P._this)

cdef extern from "<networkit/community/CoverHubDominance.hpp>":

	cdef cppclass _CoverHubDominance "NetworKit::CoverHubDominance"(_LocalCoverEvaluation):
		_CoverHubDominance(_Graph G, _Cover C) except +

cdef class CoverHubDominance(LocalCoverEvaluation):
	"""
	A quality measure that measures the dominance of hubs in clusters. The hub dominance of a single
	cluster is defined as the maximum cluster-internal degree of a node in that cluster divided by
	the maximum cluster-internal degree, i.e. the number of nodes in the cluster minus one. The
	value for all clusters is defined as the average of all clusters.
	This implementation is a natural generalization of this measure for covers.
	Strictly speaking this is not a quality measure as this is rather dependent on the type of the
	considered graph, for more information see
	Lancichinetti A, Kivel M, Saramki J, Fortunato S (2010)
	Characterizing the Community Structure of Complex Networks
	PLoS ONE 5(8): e11976. doi: 10.1371/journal.pone.0011976
	http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011976

	Parameters
	----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	C : networkit.Cover
		The cover that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _CoverHubDominance(self._G._this, self._C._this)

cdef extern from "<networkit/community/PartitionHubDominance.hpp>":

	cdef cppclass _PartitionHubDominance "NetworKit::PartitionHubDominance"(_LocalPartitionEvaluation):
		_PartitionHubDominance(_Graph G, _Partition C) except +

cdef class PartitionHubDominance(LocalPartitionEvaluation):
	"""
	A quality measure that measures the dominance of hubs in clusters. The hub dominance of a single
	cluster is defined as the maximum cluster-internal degree of a node in that cluster divided by
	the maximum cluster-internal degree, i.e. the number of nodes in the cluster minus one. The
	value for all clusters is defined as the average of all clusters.
	Strictly speaking this is not a quality measure as this is rather dependent on the type of the
	considered graph, for more information see
	Lancichinetti A, Kivel M, Saramki J, Fortunato S (2010)
	Characterizing the Community Structure of Complex Networks
	PLoS ONE 5(8): e11976. doi: 10.1371/journal.pone.0011976
	http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011976

	Parameters
	----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _PartitionHubDominance(self._G._this, self._P._this)

cdef extern from "<networkit/community/PartitionFragmentation.hpp>":

	cdef cppclass _PartitionFragmentation "NetworKit::PartitionFragmentation"(_LocalPartitionEvaluation):
		_PartitionFragmentation(_Graph G, _Partition C) except +

cdef class PartitionFragmentation(LocalPartitionEvaluation):
	"""
	This measure evaluates how fragmented a partition is. The fragmentation of a single cluster is defined as one minus the
	number of nodes in its maximum connected componented divided by its total number of nodes. Smaller values thus indicate a smaller fragmentation.

	Parameters
	----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _PartitionFragmentation(self._G._this, self._P._this)

cdef extern from "<networkit/community/StablePartitionNodes.hpp>":

	cdef cppclass _StablePartitionNodes "NetworKit::StablePartitionNodes"(_LocalPartitionEvaluation):
		_StablePartitionNodes(_Graph G, _Partition C) except +
		bool_t isStable(node u) except +

cdef class StablePartitionNodes(LocalPartitionEvaluation):
	"""
	Evaluates how stable a given partition is. A node is considered to be stable if it has strictly more connections
	to its own partition than to other partitions. Isolated nodes are considered to be stable.
	The value of a cluster is the percentage of stable nodes in the cluster.
	Larger values indicate that a clustering is more stable and thus better defined.

	Parameters
	----------
	G : networkit.Graph
		The graph on which the measure shall be evaluated
	P : networkit.Partition
		The partition that shall be evaluated
	"""
	def __cinit__(self):
		self._this = new _StablePartitionNodes(self._G._this, self._P._this)


	def isStable(self, node u):
		"""
		Check if a given node is stable, i.e. more connected to its own partition than to other partitions.

		Parameters
		----------
		u : node
			The node to check

		Returns
		-------
		bool
			If the node u is stable.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_StablePartitionNodes*>(self._this)).isStable(u)


cdef extern from "<networkit/community/CoverF1Similarity.hpp>":

	cdef cppclass _CoverF1Similarity "NetworKit::CoverF1Similarity"(_LocalCoverEvaluation):
		_CoverF1Similarity(_Graph G, _Cover C, _Cover reference) except +

cdef class CoverF1Similarity(LocalCoverEvaluation):
	"""
	Compare a given cover to a reference cover using the F1 measure.
	This is a typical similarity measure used to compare the found
	overlapping community structure to a ground truth community
	structure. Each cluster is compared to the best-matching reference
	cluster (in terms of highest F1 score). A value of 1 indicates
	perfect agreement while a while of 0 indicates complete
	disagreement. An example where this measure is used is the
	following paper:

	Alessandro Epasto, Silvio Lattanzi, and Renato Paes
	Leme. 2017. Ego-Splitting Framework: from Non-Overlapping to
	Overlapping Clusters. In Proceedings of the 23rd ACM SIGKDD
	International Conference on Knowledge Discovery and Data Mining
	(KDD '17). ACM, New York, NY, USA, 145-154. DOI:
	https://doi.org/10.1145/3097983.3098054

	Parameters
	----------
	G : Graph
		The graph on which the evaluation is performed.
	C : Cover
		The cover that shall be evaluated
        reference : Cover
		The cover to which the similarity shall be computed
	"""
	cdef Cover _reference
	def __cinit__(self, Graph G not None, Cover C not None, Cover reference not None):
		self._this = new _CoverF1Similarity(G._this, C._this, reference._this)
		self._reference = reference
		assert(self._G == G)
		assert(self._C == C)

# Module: flows

cdef extern from "<networkit/flow/EdmondsKarp.hpp>":

	cdef cppclass _EdmondsKarp "NetworKit::EdmondsKarp":
		_EdmondsKarp(const _Graph &graph, node source, node sink) except +
		void run() nogil except +
		edgeweight getMaxFlow() const
		vector[node] getSourceSet() except +
		edgeweight getFlow(node u, node v) except +
		edgeweight getFlow(edgeid eid) const
		vector[edgeweight] getFlowVector() except +

cdef class EdmondsKarp:
	"""
	The EdmondsKarp class implements the maximum flow algorithm by Edmonds and Karp.

	Parameters
	----------
	graph : networkit.Graph
		The graph
	source : node
		The source node for the flow calculation
	sink : node
		The sink node for the flow calculation
	"""
	cdef _EdmondsKarp* _this
	cdef Graph _graph

	def __cinit__(self, Graph graph not None, node source, node sink):
		self._graph = graph # store reference of graph for memory management, so the graph is not deallocated before this object
		self._this = new _EdmondsKarp(graph._this, source, sink)

	def __dealloc__(self):
		del self._this

	def run(self):
		"""
		Computes the maximum flow, executes the EdmondsKarp algorithm
		"""
		with nogil:
			self._this.run()
		return self

	def getMaxFlow(self):
		"""
		Returns the value of the maximum flow from source to sink.

		Returns
		-------
		edgeweight
			The maximum flow value
		"""
		return self._this.getMaxFlow()

	def getSourceSet(self):
		"""
		Returns the set of the nodes on the source side of the flow/minimum cut.

		Returns
		-------
		list
			The set of nodes that form the (smallest) source side of the flow/minimum cut.
		"""
		return self._this.getSourceSet()

	def getFlow(self, node u, node v = none):
		"""
		Get the flow value between two nodes u and v or an edge identified by the edge id u.
		Warning: The variant with two edge ids is linear in the degree of u.

		Parameters
		----------
		u : node or edgeid
			The first node incident to the edge or the edge id
		v : node
			The second node incident to the edge (optional if edge id is specified)

		Returns
		-------
		edgeweight
			The flow on the specified edge
		"""
		if v == none: # Assume that node and edge ids are the same type
			return self._this.getFlow(u)
		else:
			return self._this.getFlow(u, v)

	def getFlowVector(self):
		"""
		Return a copy of the flow values of all edges.

		Returns
		-------
		list
			The flow values of all edges indexed by edge id
		"""
		return self._this.getFlowVector()

# Module: properties

cdef extern from "<networkit/global/ClusteringCoefficient.hpp>" namespace "NetworKit::ClusteringCoefficient":

		double avgLocal(_Graph G, bool_t turbo) nogil except +
		double sequentialAvgLocal(_Graph G) nogil except +
		double approxAvgLocal(_Graph G, count trials) nogil except +
		double exactGlobal(_Graph G) nogil except +
		double approxGlobal(_Graph G, count trials) nogil except +

cdef class ClusteringCoefficient:
	@staticmethod
	def avgLocal(Graph G, bool_t turbo = False):
		"""
		DEPRECATED: Use centrality.LocalClusteringCoefficient and take average.

		This calculates the average local clustering coefficient of graph `G`. The graph may not contain self-loops.

		Parameters
		----------
		G : networkit.Graph
			The graph.

		Notes
		-----

		.. math:: c(G) := \\frac{1}{n} \sum_{u \in V} c(u)

		where

		.. math:: c(u) := \\frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}

		"""
		cdef double ret
		with nogil:
			ret = avgLocal(G._this, turbo)
		return ret

	@staticmethod
	def sequentialAvgLocal(Graph G):
		""" This calculates the average local clustering coefficient of graph `G` using inherently sequential triangle counting.
		Parameters
		----------
		G : networkit.Graph
			The graph.

		Notes
		-----

		.. math:: c(G) := \\frac{1}{n} \sum_{u \in V} c(u)

		where

		.. math:: c(u) := \\frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}

		"""
		cdef double ret
		with nogil:
			ret = sequentialAvgLocal(G._this)
		return ret

	@staticmethod
	def approxAvgLocal(Graph G, count trials):
		cdef double ret
		with nogil:
			ret = approxAvgLocal(G._this, trials)
		return ret

	@staticmethod
	def exactGlobal(Graph G):
		""" This calculates the global clustering coefficient. """
		cdef double ret
		with nogil:
			ret = exactGlobal(G._this)
		return ret

	@staticmethod
	def approxGlobal(Graph G, count trials):
		cdef double ret
		with nogil:
			ret = approxGlobal(G._this, trials)
		return ret

cdef extern from "<networkit/correlation/Assortativity.hpp>":

	cdef cppclass _Assortativity "NetworKit::Assortativity"(_Algorithm):
		_Assortativity(_Graph, vector[double]) except +
		_Assortativity(_Graph, _Partition) except +
		double getCoefficient() except +

cdef class Assortativity(Algorithm):
	""" """
	cdef Graph G
	cdef vector[double] attribute
	cdef Partition partition

	def __cinit__(self, Graph G, data):
		if isinstance(data, Partition):
			self._this = new _Assortativity(G._this, (<Partition>data)._this)
			self.partition = <Partition>data
		else:
			self.attribute = <vector[double]?>data
			self._this = new _Assortativity(G._this, self.attribute)
		self.G = G

	def getCoefficient(self):
		return (<_Assortativity*>(self._this)).getCoefficient()

cdef extern from "<networkit/dynamics/GraphDifference.hpp>":

	cdef cppclass _GraphDifference "NetworKit::GraphDifference"(_Algorithm):
		_GraphDifference(const _Graph &G1, const _Graph &G2) except +
		vector[_GraphEvent] getEdits() except +
		count getNumberOfEdits() except +
		count getNumberOfNodeAdditions() except +
		count getNumberOfNodeRemovals() except +
		count getNumberOfNodeRestorations() except +
		count getNumberOfEdgeAdditions() except +
		count getNumberOfEdgeRemovals() except +
		count getNumberOfEdgeWeightUpdates() except +

cdef class GraphDifference(Algorithm):
	"""
	Calculate the edge difference between two graphs.

	This calculates which graph edge additions or edge removals are
	necessary to transform one given graph into another given graph.

	Both graphs need to have the same node set, directed graphs are not
	supported currently.

	Note that edge weight differences are not detected but edge
	addition events set the correct edge weight.

	Parameters
	----------
	G1 : networkit.Graph
		The first graph to compare
	G2 : networkit.Graph
		The second graph to compare
	"""
	cdef Graph _G1, _G2

	def __cinit__(self, Graph G1, Graph G2):
		self._this = new _GraphDifference(G1._this, G2._this)
		self._G1 = G1
		self._G2 = G2

	def getEdits(self):
		""" Get the required edits.

		Returns
		-------
		list
			A list of graph events
		"""
		cdef _GraphEvent ev
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in (<_GraphDifference*>(self._this)).getEdits()]

	def getNumberOfEdits(self):
		""" Get the required number of edits.

		Returns
		-------
		int
			The number of edits.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfEdits()

	def getNumberOfNodeAdditions(self):
		""" Get the required number of node additions.

		Returns
		-------
		int
			The number of node additions.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfNodeAdditions()

	def getNumberOfNodeRemovals(self):
		""" Get the required number of node removals.

		Returns
		-------
		int
			The number of node removals.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfNodeRemovals()

	def getNumberOfNodeRestorations(self):
		""" Get the required number of node restorations.

		Returns
		-------
		int
			The number of node restorations.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfNodeRestorations()

	def getNumberOfEdgeAdditions(self):
		""" Get the required number of edge additions.

		Returns
		-------
		int
			The number of edge additions.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfEdgeAdditions()

	def getNumberOfEdgeRemovals(self):
		""" Get the required number of edge removals.

		Returns
		-------
		int
			The number of edge removals.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfEdgeRemovals()

	def getNumberOfEdgeWeightUpdates(self):
		""" Get the required number of edge weight updates.

		Returns
		-------
		int
			The number of edge weight updates.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfEdgeWeightUpdates()



# Module: dynamic

cdef extern from "<networkit/dynamics/DGSStreamParser.hpp>":

	cdef cppclass _DGSStreamParser "NetworKit::DGSStreamParser":
		_DGSStreamParser(string path, bool_t mapped, node baseIndex) except +
		vector[_GraphEvent] getStream() except +

cdef class DGSStreamParser:
	cdef _DGSStreamParser* _this

	def __cinit__(self, path, mapped=True, baseIndex=0):
		self._this = new _DGSStreamParser(stdstring(path), mapped, baseIndex)

	def __dealloc__(self):
		del self._this

	def getStream(self):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.getStream()]


cdef extern from "<networkit/dynamics/DGSWriter.hpp>":

	cdef cppclass _DGSWriter "NetworKit::DGSWriter":
		void write(vector[_GraphEvent] stream, string path) except +


cdef class DGSWriter:
	cdef _DGSWriter* _this

	def __cinit__(self):
		self._this = new _DGSWriter()

	def __dealloc__(self):
		del self._this

	def write(self, stream, path):
		cdef vector[_GraphEvent] _stream
		for ev in stream:
			_stream.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.write(_stream, stdstring(path))


# cdef extern from "<networkit/dcd2/DynamicCommunityDetection.hpp>":

# 	cdef cppclass _DynamicCommunityDetection "NetworKit::DynamicCommunityDetection":
# 		_DynamicCommunityDetection(string inputPath, string algoName, string updateStrategy, count interval, count restart, vector[string] recordSettings) except +
# 		void run() except +
# 		vector[double] getTimeline(string key) except +
# 		vector[pair[count, count]] getGraphSizeTimeline() except +
# 		vector[pair[_Graph, _Partition]] getResultTimeline() except +

# cdef class DynamicCommunityDetection:
# 	cdef _DynamicCommunityDetection* _this

# 	def __cinit__(self, inputPath, algoName, updateStrategy, interval, restart, recordSettings):
# 		self._this = new _DynamicCommunityDetection(stdstring(inputPath), stdstring(algoName), stdstring(updateStrategy), interval, restart, [stdstring(key) for key in recordSettings])

# 	def run(self):
# 		self._this.run()

# 	def getTimeline(self, key):
# 		return self._this.getTimeline(stdstring(key))

# 	def getGraphSizeTimeline(self):
# 		return self._this.getGraphSizeTimeline()

# 	def getResultTimeline(self):
# 		timeline = []
# 		for pair in self._this.getResultTimeline():
# 			_G = pair.first
# 			_zeta = pair.second
# 			timeline.append((Graph().setThis(_G), Partition().setThis(_zeta)))
# 		return timeline



cdef extern from "<networkit/dynamics/GraphUpdater.hpp>":

	cdef cppclass _GraphUpdater "NetworKit::GraphUpdater":
		_GraphUpdater(_Graph G) except +
		void update(vector[_GraphEvent] stream) nogil except +
		vector[pair[count, count]] getSizeTimeline() except +

cdef class GraphUpdater:
	""" Updates a graph according to a stream of graph events.

	Parameters
	----------
	G : networkit.Graph
	 	initial graph
	"""
	cdef _GraphUpdater* _this
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _GraphUpdater(G._this)

	def __dealloc__(self):
		del self._this

	def update(self, stream):
		cdef vector[_GraphEvent] _stream
		for ev in stream:
			_stream.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		with nogil:
			self._this.update(_stream)


# Module: coarsening

cdef extern from "<networkit/coarsening/GraphCoarsening.hpp>":

	cdef cppclass _GraphCoarsening "NetworKit::GraphCoarsening"(_Algorithm):
		_GraphCoarsening(_Graph) except +
		_Graph getCoarseGraph() except +
		vector[node] getFineToCoarseNodeMapping() except +
		map[node, vector[node]] getCoarseToFineNodeMapping() except +

cdef class GraphCoarsening(Algorithm):
	cdef Graph _G

	def __init__(self, *args, **namedargs):
		if type(self) == GraphCoarsening:
			raise RuntimeError("Error, you may not use GraphCoarsening directly, use a sub-class instead")

	def getCoarseGraph(self):
		return Graph(0).setThis((<_GraphCoarsening*>(self._this)).getCoarseGraph())

	def getFineToCoarseNodeMapping(self):
		return (<_GraphCoarsening*>(self._this)).getFineToCoarseNodeMapping()

	def getCoarseToFineNodeMapping(self):
		return (<_GraphCoarsening*>(self._this)).getCoarseToFineNodeMapping()


cdef extern from "<networkit/coarsening/ParallelPartitionCoarsening.hpp>":

	cdef cppclass _ParallelPartitionCoarsening "NetworKit::ParallelPartitionCoarsening"(_GraphCoarsening):
		_ParallelPartitionCoarsening(_Graph, _Partition, bool_t) except +


cdef class ParallelPartitionCoarsening(GraphCoarsening):
	def __cinit__(self, Graph G not None, Partition zeta not None, useGraphBuilder = True):
		self._this = new _ParallelPartitionCoarsening(G._this, zeta._this, useGraphBuilder)

cdef extern from "<networkit/coarsening/MatchingCoarsening.hpp>":

	cdef cppclass _MatchingCoarsening "NetworKit::MatchingCoarsening"(_GraphCoarsening):
		_MatchingCoarsening(_Graph, _Matching, bool_t) except +


cdef class MatchingCoarsening(GraphCoarsening):
	"""Coarsens graph according to a matching.
 	Parameters
 	----------
 	G : networkit.Graph
	M : Matching
 	noSelfLoops : bool, optional
		if true, self-loops are not produced
	"""

	def __cinit__(self, Graph G not None, Matching M not None, bool_t noSelfLoops=False):
		self._this = new _MatchingCoarsening(G._this, M._this, noSelfLoops)


# Module: scd

cdef extern from "<networkit/scd/PageRankNibble.hpp>":

	cdef cppclass _PageRankNibble "NetworKit::PageRankNibble":
		_PageRankNibble(_Graph G, double alpha, double epsilon) except +
		map[node, set[node]] run(set[node] seeds) except +

cdef class PageRankNibble:
	"""
	Produces a cut around a given seed node using the PageRank-Nibble algorithm.
	see Andersen, Chung, Lang: Local Graph Partitioning using PageRank Vectors

	Parameters:
	-----------
	G : networkit.Graph in which the cut is to be produced, must be unweighted.
	alpha : Loop probability of random walk; smaller values tend to produce larger communities.
	epsilon: Tolerance threshold for approximation of PageRank vectors
	"""
	cdef _PageRankNibble *_this
	cdef Graph _G

	def __cinit__(self, Graph G, double alpha, double epsilon):
		self._G = G
		self._this = new _PageRankNibble(G._this, alpha, epsilon)

	def run(self, set[node] seeds):
		"""
		Produces a cut around a given seed node.

		Parameters:
		-----------
		seeds : the seed node ids.
		"""
		return self._this.run(seeds)

cdef extern from "<networkit/scd/GCE.hpp>":

	cdef cppclass _GCE "NetworKit::GCE":
		_GCE(_Graph G, string quality) except +
		map[node, set[node]] run(set[node] seeds) except +

cdef class GCE:
	"""
	Produces a cut around a given seed node using the GCE algorithm.

	Parameters:
	-----------
	G : networkit.Graph in which the cut is to be produced, must be unweighted.
	"""
	cdef _GCE *_this
	cdef Graph _G

	def __cinit__(self, Graph G, quality):
		self._G = G
		self._this = new _GCE(G._this, stdstring(quality))

	def run(self, set[node] seeds):
		"""
		Produces a cut around a given seed node.

		Parameters:
		-----------
		seeds : the seed node ids.
		"""
		return self._this.run(seeds)


# Module: clique

cdef cppclass NodeVectorCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	# This is called within the run() method which is nogil!
	void cython_call_operator(const vector[node]& nodes) nogil:
		cdef bool_t error = False
		cdef string message
		# Acquire gil to allow Python code!
		with gil:
			try:
				(<object>callback)(nodes)
			except Exception as e:
				error = True
				message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
			if (error):
				throw_runtime_error(message)

cdef extern from "<networkit/clique/MaximalCliques.hpp>":

	cdef cppclass _MaximalCliques "NetworKit::MaximalCliques"(_Algorithm):
		_MaximalCliques(_Graph G, bool_t maximumOnly) except +
		_MaximalCliques(_Graph G, NodeVectorCallbackWrapper callback) except +
		vector[vector[node]] getCliques() except +

cdef class MaximalCliques(Algorithm):
	"""
	Algorithm for listing all maximal cliques.

	The implementation is based on the "hybrid" algorithm described in

	Eppstein, D., & Strash, D. (2011).
	Listing All Maximal Cliques in Large Sparse Real-World Graphs.
	In P. M. Pardalos & S. Rebennack (Eds.),
	Experimental Algorithms (pp. 364375). Springer Berlin Heidelberg.
	Retrieved from http://link.springer.com/chapter/10.1007/978-3-642-20662-7_31

	The running time of this algorithm should be in O(d^2 * n * 3^{d/3})
	where f is the degeneracy of the graph, i.e., the maximum core number.
	The running time in practive depends on the structure of the graph. In
	particular for complex networks it is usually quite fast, even graphs with
	millions of edges can usually be processed in less than a minute.

	Parameters
	----------
	G : networkit.Graph
		The graph to list the cliques for
	maximumOnly : bool
		A value of True denotes that only one maximum clique is desired. This enables
		further optimizations of the algorithm to skip smaller cliques more
		efficiently. This parameter is only considered when no callback is given.
	callback : callable
		If a callable Python object is given, it will be called once for each
		maximal clique. Then no cliques will be stored. The callback must accept
		one parameter which is a list of nodes.
	"""
	cdef NodeVectorCallbackWrapper* _callback
	cdef Graph _G
	cdef object _py_callback

	def __cinit__(self, Graph G not None, bool_t maximumOnly = False, object callback = None):
		self._G = G

		if callable(callback):
			# Make sure the callback is not de-allocated!
			self._py_callback = callback
			self._callback = new NodeVectorCallbackWrapper(callback)
			try:
				self._this = new _MaximalCliques(self._G._this, dereference(self._callback))
			except BaseException as e:
				del self._callback
				self._callback = NULL
				raise e
		else:
			self._callback = NULL
			self._this = new _MaximalCliques(self._G._this, maximumOnly)

	def __dealloc__(self):
		if not self._callback == NULL:
			del self._callback
			self._callback = NULL

	def getCliques(self):
		"""
		Return all found cliques unless a callback was given.

		This method will throw if a callback was given and thus the cliques were not stored.
		If only the maximum clique was stored, it will return exactly one clique unless the graph
		is empty.

		Returns
		-------
		A list of cliques, each being represented as a list of nodes.
		"""
		return (<_MaximalCliques*>(self._this)).getCliques()

# Module: linkprediction

cdef extern from "<networkit/linkprediction/LinkPredictor.hpp>":

	cdef cppclass _LinkPredictor "NetworKit::LinkPredictor":
		_LinkPredictor(const _Graph& G) except +
		double run(node u, node v) except +
		vector[pair[pair[node, node], double]] runAll() except +
		vector[pair[pair[node, node], double]] runOn(vector[pair[node, node]] nodePairs) except +
		void setGraph(const _Graph& newGraph) except +

cdef class LinkPredictor:
	""" Abstract base class for link predictors.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""
	cdef _LinkPredictor* _this

	def __cinit__(self, *args):
		# The construction is handled by the subclasses
		return

	def __dealloc__(self):
		if self._this is not NULL:
			del self._this
			self._this = NULL

	def setGraph(self, Graph newGraph):
		""" Sets the graph to work on.

		Parameters
		----------
		newGraph : networkit.Graph
			The graph to work on.
   	"""
		self._this.setGraph(newGraph._this)

	def run(self, node u, node v):
		""" Returns a score indicating the likelihood of a future link between the given nodes.

		Prior to calling this method a graph should be provided through the constructor or
		by calling setGraph. Note that only undirected graphs are accepted.
		There is also no lower or upper bound for scores and the actual range of values depends
		on the specific link predictor implementation. In case u == v a 0 is returned.
		If suitable this method might make use of parallelization to enhance performance.

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		A prediction-score indicating the likelihood of a future link between the given nodes.
		"""
		return self._this.run(u, v)

	def runAll(self):
		""" Runs the link predictor on all currently unconnected node-pairs.

		Possible self-loops are also excluded. The method makes use of parallelisation.

		Returns
		-------
		A vector of pairs containing all currently unconnected node-pairs as the first elements
		and the corresponding scores as the second elements. The vector is sorted ascendingly by node-pair.
		"""
		return move(self._this.runAll())

	def runOn(self, vector[pair[node, node]] nodePairs):
		""" Executes the run-method on al given node-pairs and returns a vector of predictions.

		The result is a vector of pairs where the first element is the node-pair and it's second
		element the corresponding score generated by the run-method. The method makes use of
		parallelisation.

		Parameters
		----------
		nodePairs : vector[pair[node, node]]
			Node-pairs to run the predictor on.

		Returns
		-------
		A vector of pairs containing the given node-pair as the first element and it's corresponding score
		as the second element. The vector is sorted ascendingly by node-pair.
		"""
		return move(self._this.runOn(nodePairs))

cdef extern from "<networkit/linkprediction/KatzIndex.hpp>":

	cdef cppclass _KatzIndex "NetworKit::KatzIndex"(_LinkPredictor):
		_KatzIndex(count maxPathLength, double dampingValue) except +
		_KatzIndex(const _Graph& G, count maxPathLength, double dampingValue) except +

cdef class KatzIndex(LinkPredictor):
	""" Implementation of the Katz index.

	Katz index assigns a pair of nodes a similarity score
	that is based on the sum of the weighted number of paths of length l
	where l is smaller than a given limit.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to operate on. Defaults to None.
	maxPathLength : count, optional
		Maximal length of the paths to consider. Defaults to 5.
	dampingValue : double, optional
		Used to exponentially damp every addend of the sum. Should be in (0, 1]. Defaults to 0.005.
	"""

	def __cinit__(self, Graph G = None, count maxPathLength = 5, double dampingValue = 0.005):
		if G is None:
			self._this = new _KatzIndex(maxPathLength, dampingValue)
		else:
			self._this = new _KatzIndex(G._this, maxPathLength, dampingValue)

	def run(self, node u, node v):
		""" Returns the similarity score for the given node-pair based on the Katz index specified during construction.

		The algorithm considers all paths starting at the node with the smaller degree except the algorithm
		started at the other node at the last call.

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The similarity score of the given node-pair calculated by the specified Katz index.
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/CommonNeighborsIndex.hpp>":

	cdef cppclass _CommonNeighborsIndex "NetworKit::CommonNeighborsIndex"(_LinkPredictor):
		_CommonNeighborsIndex() except +
		_CommonNeighborsIndex(const _Graph& G) except +

cdef class CommonNeighborsIndex(LinkPredictor):
	""" The CommonNeighborsIndex calculates the number of common neighbors of a node-pair in a given graph.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _CommonNeighborsIndex()
		else:
			self._this = new _CommonNeighborsIndex(G._this)

	def run(self, node u, node v):
		""" Returns the number of common neighbors of the given nodes u and v.

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The number of common neighbors of u and v.
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/PreferentialAttachmentIndex.hpp>":

	cdef cppclass _PreferentialAttachmentIndex "NetworKit::PreferentialAttachmentIndex"(_LinkPredictor):
		_PreferentialAttachmentIndex() except +
		_PreferentialAttachmentIndex(const _Graph& G) except +

cdef class PreferentialAttachmentIndex(LinkPredictor):
	""" Implementation of the Preferential Attachment Index.

	The run-method simply calculates the product of the number of nodes in the neighborhoods
	regarding the given nodes.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _PreferentialAttachmentIndex()
		else:
			self._this = new _PreferentialAttachmentIndex(G._this)

	def run(self, node u, node v):
		""" Returns the product of the cardinalities of the neighborhoods regarding u and v.

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The product of the cardinalities of the neighborhoods regarding u and v
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/JaccardIndex.hpp>":

	cdef cppclass _JaccardIndex "NetworKit::JaccardIndex"(_LinkPredictor):
		_JaccardIndex() except +
		_JaccardIndex(const _Graph& G) except +

cdef class JaccardIndex(LinkPredictor):
	""" Implementation of the Jaccard index which normalizes the Common Neighbors Index.

	This is done through dividing the number of common neighbors by the number of nodes
	in the neighboorhood-union.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""
	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _JaccardIndex()
		else:
			self._this = new _JaccardIndex(G._this)

	def run(self, node u, node v):
		""" Returns the Jaccard index for the given node-pair (u, v).

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The Jaccard index for the given node-pair (u, v).
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/AdamicAdarIndex.hpp>":

	cdef cppclass _AdamicAdarIndex "NetworKit::AdamicAdarIndex"(_LinkPredictor):
		_AdamicAdarIndex() except +
		_AdamicAdarIndex(const _Graph& G) except +

cdef class AdamicAdarIndex(LinkPredictor):
	""" Implementation of the Adamic/Adar Index.

	The index sums up the reciprocals of the logarithm of the degree of all
	common neighbors of u and v.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _AdamicAdarIndex()
		else:
			self._this = new _AdamicAdarIndex(G._this)

	def run(self, node u, node v):
		""" Returns the Adamic/Adar Index of the given node-pair (u, v).

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The Adamic/Adar Index of the given node-pair (u, v).
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/UDegreeIndex.hpp>":

	cdef cppclass _UDegreeIndex "NetworKit::UDegreeIndex"(_LinkPredictor):
		_UDegreeIndex() except +
		_UDegreeIndex(const _Graph& G) except +

cdef class UDegreeIndex(LinkPredictor):
	""" Index that simply returns the degree of the first given node.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _UDegreeIndex()
		else:
			self._this = new _UDegreeIndex(G._this)

	def run(self, node u, node v):
		""" Returns the degree of the first node provided, namely u.

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The degree of the first node provided, namely u.
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/VDegreeIndex.hpp>":

	cdef cppclass _VDegreeIndex "NetworKit::VDegreeIndex"(_LinkPredictor):
		_VDegreeIndex() except +
		_VDegreeIndex(const _Graph& G) except +

cdef class VDegreeIndex(LinkPredictor):
	""" Index that simply returns the degree of the second given node.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _VDegreeIndex()
		else:
			self._this = new _VDegreeIndex(G._this)

	def run(self, node u, node v):
		""" Returns the degree of the second node provided, namely v.

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The degree of the second node provided, namely v.
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/AlgebraicDistanceIndex.hpp>":

	cdef cppclass _AlgebraicDistanceIndex "NetworKit::AlgebraicDistanceIndex"(_LinkPredictor):
		_AlgebraicDistanceIndex(count numberSystems, count numberIterations, double omega, index norm) except +
		_AlgebraicDistanceIndex(const _Graph& G, count numberSystems, count numberIterations, double omega, index norm) except +
		void preprocess() except +
		double run(node u, node v) except +

cdef class AlgebraicDistanceIndex(LinkPredictor):
	""" Algebraic distance assigns a distance value to pairs of nodes according to their structural closeness in the graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to work on. Can be set to None and default is None.
	numberSystems : count
		Number of vectors/systems used for algebraic iteration.
	numberIterations : count
		Number of iterations in each system.
	omega : double, optional
		Overrelaxation parameter, default: 0.5.
	norm : index, optional
		The norm factor of the extended algebraic distance. Maximum norm is realized by setting the norm to 0. Default: 2.
	"""

	def __cinit__(self, Graph G, count numberSystems, count numberIterations, double omega = 0.5, index norm = 2):
		if G is None:
			self._this = new _AlgebraicDistanceIndex(numberSystems, numberIterations, omega, norm)
		else:
			self._this = new _AlgebraicDistanceIndex(G._this, numberSystems, numberIterations, omega, norm)

	def preprocess(self):
		""" Executes necessary initializations.

		Starting with random initialization, compute for all numberSystems
		"diffusion" systems the situation after numberIterations iterations
		of overrelaxation with overrelaxation parameter omega.

		REQ: Needs to be called before algdist delivers meaningful results!
		"""
		(<_AlgebraicDistanceIndex *>self._this).preprocess()

	def run(self, node u, node v):
		""" Returns the extended algebraic distance between node u and node v in the norm specified in the constructor.

		Parameters
		----------
		u : node
			The first node.
		v : node
			The second node.

		Returns
		-------
		Extended algebraic distance between the two nodes.
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/NeighborhoodDistanceIndex.hpp>":

	cdef cppclass _NeighborhoodDistanceIndex "NetworKit::NeighborhoodDistanceIndex"(_LinkPredictor):
		_NeighborhoodDistanceIndex() except +
		_NeighborhoodDistanceIndex(const _Graph& G) except +
		double run(node u, node v) except +

cdef class NeighborhoodDistanceIndex(LinkPredictor):
	""" Assigns a distance value to pairs of nodes according to the overlap of their neighborhoods.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""
	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _NeighborhoodDistanceIndex()
		else:
			self._this = new _NeighborhoodDistanceIndex(G._this)

	def run(self, node u, node v):
		""" Returns the Neighborhood Distance index for the given node-pair (u, v).

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The Neighborhood Distance index for the given node-pair (u, v).
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/TotalNeighborsIndex.hpp>":

	cdef cppclass _TotalNeighborsIndex "NetworKit::TotalNeighborsIndex"(_LinkPredictor):
		_TotalNeighborsIndex() except +
		_TotalNeighborsIndex(const _Graph& G) except +

cdef class TotalNeighborsIndex(LinkPredictor):
	""" Implementation of the Total Neighbors Index.

	This index is also known as Total Friends Index and returns
	the number of nodes in the neighborhood-union of u and v.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _TotalNeighborsIndex()
		else:
			self._this = new _TotalNeighborsIndex(G._this)

	def run(self, node u, node v):
		""" Returns the number of total union-neighbors for the given node-pair (u, v).

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The number of total union-neighbors for the given node-pair (u, v).
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/NeighborsMeasureIndex.hpp>":

	cdef cppclass _NeighborsMeasureIndex "NetworKit::NeighborsMeasureIndex"(_LinkPredictor):
		_NeighborsMeasureIndex() except +
		_NeighborsMeasureIndex(const _Graph& G) except +

cdef class NeighborsMeasureIndex(LinkPredictor):
	""" Implementation of the Neighbors Measure Index.

	This index is also known as Friends Measure and simply returns
	the number of connections between neighbors of the given nodes u and v.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _NeighborsMeasureIndex()
		else:
			self._this = new _NeighborsMeasureIndex(G._this)

	def run(self, node u, node v):
		""" Returns the number of connections between neighbors of u and v.

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The number of connections between neighbors of u and v.
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/SameCommunityIndex.hpp>":

	cdef cppclass _SameCommunityIndex "NetworKit::SameCommunityIndex"(_LinkPredictor):
		_SameCommunityIndex() except +
		_SameCommunityIndex(const _Graph& G) except +

cdef class SameCommunityIndex(LinkPredictor):
	""" Index to determine whether two nodes are in the same community.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _SameCommunityIndex()
		else:
			self._this = new _SameCommunityIndex(G._this)

	def run(self, node u, node v):
		""" Returns 1 if the given nodes u and v are in the same community, 0 otherwise.

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		1 if the given nodes u and v are in the same community, 0 otherwise.
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/AdjustedRandIndex.hpp>":

	cdef cppclass _AdjustedRandIndex "NetworKit::AdjustedRandIndex"(_LinkPredictor):
		_AdjustedRandIndex() except +
		_AdjustedRandIndex(const _Graph& G) except +

cdef class AdjustedRandIndex(LinkPredictor):
	""" AdjustedRandIndex proposed by Hoffman et al. with natural threshold of 0.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _AdjustedRandIndex()
		else:
			self._this = new _AdjustedRandIndex(G._this)

	def run(self, node u, node v):
		""" Returns the Adjusted Rand Index of the given node-pair (u, v).

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The Adjusted Rand Index of the given node-pair (u, v).
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/ResourceAllocationIndex.hpp>":

	cdef cppclass _ResourceAllocationIndex "NetworKit::ResourceAllocationIndex"(_LinkPredictor):
		_ResourceAllocationIndex() except +
		_ResourceAllocationIndex(const _Graph& G) except +

cdef class ResourceAllocationIndex(LinkPredictor):
	""" Implementation of the ResourceAllocationIndex.

	The index is similar to Adamic/Adar and sums up the reciprocals of
	the degree of all common neighbors of u and v.

	Parameters
	----------
	G : networkit.Graph, optional
		The graph to work on. Defaults to None.
	"""

	def __cinit__(self, Graph G = None):
		if G is None:
			self._this = new _ResourceAllocationIndex()
		else:
			self._this = new _ResourceAllocationIndex(G._this)

	def run(self, node u, node v):
		""" Returns the Resource Allocation Index of the given node-pair (u, v).

		Parameters
		----------
		u : node
			First node in graph.
		v : node
			Second node in graph.

		Returns
		-------
		The Resource Allocation Index of the given node-pair (u, v).
		"""
		return self._this.run(u, v)

cdef extern from "<networkit/linkprediction/RandomLinkSampler.hpp>" namespace "NetworKit::RandomLinkSampler":

	_Graph byPercentage(_Graph G, double percentage) except +
	_Graph byCount(_Graph G, count numLinks) except +

cdef class RandomLinkSampler:
	""" Provides methods to randomly sample a number of edges from a given graph. """

	@staticmethod
	def byPercentage(Graph G, double percentage):
		""" Returns a graph that contains percentage percent of links form the given graph G.

		The links are randomly selected from G until the given percentage is reached.

		Parameters
		----------
		G : networkit.Graph
			The graph to construct the training graph from.
		percentage : double
			Percentage of links regarding the number of links in the given graph that should
			be in the returned graph.

		Returns
		-------
		A graph that contains the given percentage of links from G.
		"""
		return Graph().setThis(byPercentage(G._this, percentage))

	@staticmethod
	def byCount(Graph G, count numLinks):
		""" Returns a graph that contains numLinks links from the given graph G.

		The links are randomly selected from G until the given count is reached.

		Parameters
		----------
		G : networkit.Graph
			The graph to construct the training graph from.
		numLinks : count
			Number of links the returned graph should consist of.

		Returns
		-------
		A graph that contains the given number of links from G.
		"""
		return Graph().setThis(byCount(G._this, numLinks))

cdef extern from "<networkit/linkprediction/EvaluationMetric.hpp>":

	cdef cppclass _EvaluationMetric "NetworKit::EvaluationMetric":
		_EvaluationMetric() except +
		_EvaluationMetric(const _Graph& testGraph) except +
		void setTestGraph(const _Graph& newTestGraph) except +
		pair[vector[double], vector[double]] getCurve(vector[pair[pair[node, node], double]] predictions, count numThresholds) except +
		double getAreaUnderCurve() except +
		double getAreaUnderCurve(pair[vector[double], vector[double]] curve) except +

cdef class EvaluationMetric:
	""" Abstract base class for evaluation curves.

	The evualation curves are generated based on the predictions calculated
	by the link predictor and a testGraph to compare against.

	Parameters
	----------
	testGraph : networkit.Graph
		Graph containing the links to use for evaluation. Can be set to None and default is None.
	"""
	cdef _EvaluationMetric *_this

	def __cinit__(self, *args):
		# The construction is handled by the subclasses
		return

	def __dealloc__(self):
		if self._this is not NULL:
			del self._this
			self._this = NULL

	def setTestGraph(self, Graph newTestGraph):
		""" Sets a new graph to use as ground truth for evaluation.

		Note that this won't reset the most recently calculated curve and as a consequence
		getAreaUnderCurve() will still behave as expected by returning the AUC of the most recent curve.

		Parameters
		----------
		newTestGraph : networkit.Graph
			New graph to use as ground truth.
		"""
		self._this.setTestGraph(newTestGraph._this)

	def getCurve(self, vector[pair[pair[node, node], double]] predictions, count numThresholds = 1000):
		""" Returns a pair of X- and Y-vectors describing the evaluation curve generated from the given predictions.

		The latest y-value will be used as a tie-breaker in case there are multiple y-values for one x-value.
		Note that the given number of thresholds (@a numThresholds) is an upper bound for the number of
		points returned. This is due to the fact that multiple y-values can map to one x-value in which case
		the tie-breaking behaviour described above will intervene.

		Parameters
		----------
		predictions : vector[pair[pair[node, node], double]]
			Predictions to evaluate.
		numThresholds : count, optional
			The number of thresholds to use the metric on. Defaults to 1000.

		Returns
		-------
		A pair of vectors where the first vectors contains all x-values and the second one contains the
		corresponding y-value.
		"""
		return self._this.getCurve(predictions, numThresholds)

	def getAreaUnderCurve(self, pair[vector[double], vector[double]] curve = pair[vector[double], vector[double]]()):
		""" Returns the area under the most recently calculated or optionally the given curve by using the trapezoidal rule.

		Note that if there is no curve specified or the vectors of the given curves are empty than
		the area under the most recently calculated curve will be returned.

		Parameters
		----------
		curve : pair[vector[double], vector[double]]
			Curve whose AUC to determine. Default: Pair of empty vectors.

		Returns
		-------
		The area under the given curve.
		"""
		if len(curve.first) == 0:
			return self._this.getAreaUnderCurve()
		return self._this.getAreaUnderCurve(curve)

cdef extern from "<networkit/linkprediction/ROCMetric.hpp>":

	cdef cppclass _ROCMetric "NetworKit::ROCMetric"(_EvaluationMetric):
		_ROCMetric() except +
		_ROCMetric(const _Graph& testGraph) except +
		pair[vector[double], vector[double]] getCurve(vector[pair[pair[node, node], double]] predictions, count numThresholds) except +

cdef class ROCMetric(EvaluationMetric):
	""" Provides points that define the Receiver Operating Characteristic curve for a given set of predictions.

	Based on the generated points the area under the curve can be calculated with the trapzoidal rule.

	Parameters
	----------
	testGraph : networkit.Graph, optional
		Graph containing the links to use for evaluation. Defaults to None.
	"""

	def __cinit__(self, Graph testGraph = None):
		if testGraph is None:
			self._this = new _ROCMetric()
		else:
			self._this = new _ROCMetric(testGraph._this)

	def getCurve(self, vector[pair[pair[node, node], double]] predictions, count numThresholds = 1000):
		""" Generate the points of the Receiver Operating Characteristic curve regarding the previously set predictions.

		Note that in the case of multiple y-values mapping to the same x-value the highest (=latest) y-value gets picked.

		Parameters
		----------
		predictions : vector[pair[pair[node, node], double]]
			Predictions to evaluate.
		numThresholds : count, optional
			The number of thresholds to use the metric on. Defaults to 1000.

		Returns
		-------
		A pair of vectors where the first vector contains the false positive rates and the second vector the
		corresponding true positive rates.
		"""
		return self._this.getCurve(predictions, numThresholds)

cdef extern from "<networkit/linkprediction/PrecisionRecallMetric.hpp>":

	cdef cppclass _PrecisionRecallMetric "NetworKit::PrecisionRecallMetric"(_EvaluationMetric):
		_PrecisionRecallMetric() except +
		_PrecisionRecallMetric(const _Graph& testGraph) except +
		pair[vector[double], vector[double]] getCurve(vector[pair[pair[node, node], double]] predictions, count numThresholds) except +

cdef class PrecisionRecallMetric(EvaluationMetric):
	""" Provides points that define the Precision-Recall curve for a given set of predictions.

	Based on the generated points the area under the curve can be calculated with the trapzoidal rule.

	Parameters
	----------
	testGraph : networkit.Graph, optional
		Graph containing the links to use for evaluation. Defaults to None.
	"""

	def __cinit__(self, Graph testGraph = None):
		if testGraph is None:
			self._this = new _PrecisionRecallMetric()
		else:
			self._this = new _PrecisionRecallMetric(testGraph._this)

	def getCurve(self, vector[pair[pair[node, node], double]] predictions, count numThresholds = 1000):
		""" Generates the points for the Precision-Recall curve with respect to the given predictions.

		The curve assigns every recall-value a corresponding precision as the y-value.
		In case of a tie regarding multiple y-values for a x-value the smallest (= latest) y-value will be used.

		Parameters
		----------
		predictions : vector[pair[pair[node, node], double]]
			Predictions to evaluate.
		numThresholds : count, optional
			The number of thresholds to use the metric on. Defaults to 1000.

		Returns
		-------
		A pair of vectors where the first vector contains all recall-values and the second vector
		the corresponding precision-values.
		"""
		return self._this.getCurve(predictions, numThresholds)

cdef extern from "<networkit/linkprediction/MissingLinksFinder.hpp>":

	cdef cppclass _MissingLinksFinder "NetworKit::MissingLinksFinder":
		_MissingLinksFinder(const _Graph& G) except +
		vector[pair[node, node]] findAtDistance(count k) except +
		vector[pair[node, node]] findFromNode(node u, count k) except +

cdef class MissingLinksFinder:
	""" Allows the user to find missing links in the given graph.

	The absent links to find are narrowed down by providing a distance
	that the nodes of the missing links should have.
	For example in case of distance 2 only node-pairs that would close
	a triangle in the given graph get returned.

	Parameters
	----------
	G : networkit.Graph
		The graph to find missing links in.
	"""
	cdef _MissingLinksFinder* _this

	def __cinit__(self, Graph G):
		self._this = new _MissingLinksFinder(G._this)

	def __dealloc__(self):
		del self._this

	def findAtDistance(self, count k):
		""" Returns all missing links in the graph that have distance k.

		Note that a distance of k actually means that there are k different links
		on the path of the two nodes that are connected through that path.

		Parameters
		----------
		k : count
			Distance of the absent links.

		Returns
		-------
		An ascendingly sorted vector of node-pairs where there is a missing link of distance k
		between the two nodes.
		"""
		return move(self._this.findAtDistance(k))

	def findFromNode(self, node u, count k):
		""" Returns all missing links in the graph that have distance k and are connected to u.

		Note that a distance of k actually means that there are k different links
		on the path of the two nodes that are connected through that path.

		Parameters
		----------
		u : node
			Node to find missing links from.
		k : count
			Distance of the absent links.

		Returns
		-------
		A vector of node-pairs where there is a missing link of distance k
		between the given node u and another node in the graph.
		"""
		return move(self._this.findFromNode(u, k))

cdef extern from "<networkit/linkprediction/NeighborhoodUtility.hpp>" namespace "NetworKit::NeighborhoodUtility":

	vector[node] getNeighborsUnion(const _Graph& G, node u, node v) except +
	vector[node] getCommonNeighbors(const _Graph& G, node u, node v) except +

cdef class NeighborhoodUtility:
	""" Provides basic operations on neighborhoods in a given graph. """

	@staticmethod
	def getNeighborsUnion(Graph G, node u, node v):
		""" Returns the union of the neighboorhoods of u and v.

		Parameters
		----------
		G : networkit.Graph
			Graph to obtain neighbors-union from.
		u : node
			First node.
		v : node
			Second node.

		Returns
		-------
		A vector containing all the nodes in the neighboorhood-union of u and v.
		"""
		return getNeighborsUnion(G._this, u, v)

	@staticmethod
	def getCommonNeighbors(Graph G, node u, node v):
		""" Returns a vector containing the node-ids of all common neighbors of u and v.

		Parameters
		----------
		G : networkit.Graph
			Graph to obtain common neighbors from.
		u : node
			First node.
		v : node
			Second node.

		Returns
		-------
		A vector containing the node-ids of all common neighbors of u and v.
		"""
		return getCommonNeighbors(G._this, u, v)

cdef extern from "<networkit/linkprediction/LinkThresholder.hpp>" namespace "NetworKit::LinkThresholder":

	vector[pair[node, node]] byScore(vector[pair[pair[node, node], double]] predictions, double minScore)
	vector[pair[node, node]] byCount(vector[pair[pair[node, node], double]] predictions, count numLinks)
	vector[pair[node, node]] byPercentage(vector[pair[pair[node, node], double]] predictions, double percentageLinks)

cdef class LinkThresholder:
	""" Filters given predictions based on some criterion and returns a vector of node-pairs that fulfill the given criterion.

	This can be used to determine which node-pairs should actually be interpreted
	as future links and which shouldn't.
	"""

	@staticmethod
	def byScore(vector[pair[pair[node, node], double]] predictions, double minScore):
		""" Returns the node-pairs whose scores are at least equal to the given minScore.

		Parameters
		----------
		predictions : vector[pair[pair[node, node], double]].
			Predictions to filter.
		minScore : double
			Minimal score that the returned node-pairs should have.

		Returns
		-------
		A vector of node-pairs whose scores are at least equal to the given minScore.
		"""
		return byScore(predictions, minScore)

	@staticmethod
	def byCount(vector[pair[pair[node, node], double]] predictions, count numLinks):
		""" Returns the first numLinks highest scored node-pairs.

		Parameters
		----------
		predictions : vector[pair[pair[node, node], double]].
			Predictions to filter.
		numLinks : count
			Number of top-scored node-pairs to return.

		Returns
		-------
		The first numLinks highest scored node-pairs.
		"""
		return byCount(predictions, numLinks)

	@staticmethod
	def byPercentage(vector[pair[pair[node, node], double]] predictions, double percentageLinks):
		""" Returns the first percentageLinks percent of the highest scores node-pairs.

		Parameters
		----------
		predictions : vector[pair[pair[node, node], double]].
			Predictions to filter.
		percentageLinks : double
			Percentage of highest scored node-pairs to return.

		Returns
		-------
		The first percentageLinks percent of the highest scores node-pairs.
		"""
		return byPercentage(predictions, percentageLinks)

cdef extern from "<networkit/linkprediction/PredictionsSorter.hpp>" namespace "NetworKit::PredictionsSorter":

	void sortByScore (vector[pair[pair[node, node], double]]& predictions) except +
	void sortByNodePair (vector[pair[pair[node, node], double]]& predictions) except +

cdef class PredictionsSorter:
	""" Allows the sorting of predictions by score or node-pair. """

	@staticmethod
	def sortByScore(list predictions):
		""" Sorts the given predictions descendingly by score.

		In case there is a tie the node-pairs are used as a tie-breaker by sorting them
		ascendingly on the first node and on a tie ascendingly by the second node.

		Parameters
		----------
		predictions : vector[pair[pair[node, node], double]]
			The predictions to sort.
		"""
		cdef vector[pair[pair[node, node], double]] predCopy = predictions
		sortByScore(predCopy)
		predictions[:] = predCopy

	@staticmethod
	def sortByNodePair(list predictions):
		""" Sorts the predictions ascendingly by node-pair.

		This means for example (0, 0) < (0, 1) and (1, 1) < (1, 0).

		Parameters
		----------
		predictions : vector[pair[pair[node, node], double]]
			The predictions to sort.
		"""
		cdef vector[pair[pair[node, node], double]] predCopy = predictions
		sortByNodePair(predCopy)
		predictions[:] = predCopy

# Module: EdgeScore

cdef extern from "<networkit/edgescores/EdgeScore.hpp>":

	cdef cppclass _EdgeScore "NetworKit::EdgeScore"[T](_Algorithm):
		_EdgeScore(const _Graph& G) except +
		vector[T] scores() except +
		T score(edgeid eid) except +
		T score(node u, node v) except +

cdef class EdgeScore(Algorithm):
	"""
	TODO DOCSTIRNG
	"""
	cdef Graph _G

	cdef bool_t isDoubleValue(self):
		raise RuntimeError("Implement in subclass")

	def __init__(self, *args, **namedargs):
		if type(self) == EdgeScore:
			raise RuntimeError("Error, you may not use EdgeScore directly, use a sub-class instead")

	def __dealloc__(self):
		self._G = None # just to be sure the graph is deleted

	def score(self, u, v = None):
		if v is None:
			if self.isDoubleValue():
				return (<_EdgeScore[double]*>(self._this)).score(u)
			else:
				return (<_EdgeScore[count]*>(self._this)).score(u)
		else:
			if self.isDoubleValue():
				return (<_EdgeScore[double]*>(self._this)).score(u, v)
			else:
				return (<_EdgeScore[count]*>(self._this)).score(u, v)

	def scores(self):
		if self.isDoubleValue():
			return (<_EdgeScore[double]*>(self._this)).scores()
		else:
			return (<_EdgeScore[count]*>(self._this)).scores()


cdef extern from "<networkit/edgescores/ChibaNishizekiTriangleEdgeScore.hpp>":

	cdef cppclass _ChibaNishizekiTriangleEdgeScore "NetworKit::ChibaNishizekiTriangleEdgeScore"(_EdgeScore[count]):
		_ChibaNishizekiTriangleEdgeScore(const _Graph& G) except +

cdef class ChibaNishizekiTriangleEdgeScore(EdgeScore):
	"""
	Calculates for each edge the number of triangles it is embedded in.

	Parameters
	----------
	G : networkit.Graph
		The graph to count triangles on.
	"""

	def __cinit__(self, Graph G):
		"""
		G : networkit.Graph
			The graph to count triangles on.
		"""
		self._G = G
		self._this = new _ChibaNishizekiTriangleEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return False

cdef extern from "<networkit/edgescores/ChibaNishizekiQuadrangleEdgeScore.hpp>":

	cdef cppclass _ChibaNishizekiQuadrangleEdgeScore "NetworKit::ChibaNishizekiQuadrangleEdgeScore"(_EdgeScore[count]):
		_ChibaNishizekiQuadrangleEdgeScore(const _Graph& G) except +

cdef class ChibaNishizekiQuadrangleEdgeScore(EdgeScore):
	"""
	Calculates for each edge the number of quadrangles (circles of length 4) it is embedded in.

	Parameters
	----------
	G : networkit.Graph
		The graph to count quadrangles on.
	"""

	def __cinit__(self, Graph G):
		"""
		Parameters
		----------
		G : networkit.Graph
			The graph to count quadrangles on.
		"""
		self._G = G
		self._this = new _ChibaNishizekiQuadrangleEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return False

cdef extern from "<networkit/edgescores/TriangleEdgeScore.hpp>":

	cdef cppclass _TriangleEdgeScore "NetworKit::TriangleEdgeScore"(_EdgeScore[double]):
		_TriangleEdgeScore(const _Graph& G) except +

cdef class TriangleEdgeScore(EdgeScore):
	"""
	Triangle counting.

	Parameters
	----------
	G : networkit.Graph
		The graph to count triangles on.
	"""

	def __cinit__(self, Graph G):
		"""
		Parameters
		----------
		G : networkit.Graph
			The graph to count triangles on.
		"""
		self._G = G
		self._this = new _TriangleEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return False

cdef extern from "<networkit/edgescores/EdgeScoreLinearizer.hpp>":

	cdef cppclass _EdgeScoreLinearizer "NetworKit::EdgeScoreLinearizer"(_EdgeScore[double]):
		_EdgeScoreLinearizer(const _Graph& G, const vector[double]& attribute, bool_t inverse) except +

cdef class EdgeScoreLinearizer(EdgeScore):
	"""
	Linearizes a score such that values are evenly distributed between 0 and 1.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	a : vector[double]
		Edge score that shall be linearized.
	"""
	cdef vector[double] _score

	def __cinit__(self, Graph G, vector[double] score, inverse = False):
		self._G = G
		self._score = score
		self._this = new _EdgeScoreLinearizer(G._this, self._score, inverse)

	cdef bool_t isDoubleValue(self):
		return True


cdef extern from "<networkit/edgescores/EdgeScoreNormalizer.hpp>":

	cdef cppclass _EdgeScoreNormalizer "NetworKit::EdgeScoreNormalizer"[T](_EdgeScore[double]):
		_EdgeScoreNormalizer(const _Graph&, vector[T]&, bool_t inverse, double lower, double upper) except +

cdef class EdgeScoreNormalizer(EdgeScore):
	"""
	Normalize an edge score such that it is in a certain range.

	Parameters
	----------
	G : networkit.Graph
		The graph the edge score is defined on.
	score : vector[double]
		The edge score to normalize.
	inverse
		Set to True in order to inverse the resulting score.
	lower
		Lower bound of the target range.
	upper
		Upper bound of the target range.
	"""
	cdef vector[double] _inScoreDouble
	cdef vector[count] _inScoreCount

	def __cinit__(self, Graph G not None, score, bool_t inverse = False, double lower = 0.0, double upper = 1.0):
		self._G = G
		try:
			self._inScoreDouble = <vector[double]?>score
			self._this = new _EdgeScoreNormalizer[double](G._this, self._inScoreDouble, inverse, lower, upper)
		except TypeError:
			try:
				self._inScoreCount = <vector[count]?>score
				self._this = new _EdgeScoreNormalizer[count](G._this, self._inScoreCount, inverse, lower, upper)
			except TypeError:
				raise TypeError("score must be either a vector of integer or float")

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/edgescores/EdgeScoreBlender.hpp>":

	cdef cppclass _EdgeScoreBlender "NetworKit::EdgeScoreBlender"(_EdgeScore[double]):
		_EdgeScoreBlender(const _Graph&, const vector[double]&, const vector[double]&, const vector[bool_t]&) except +

cdef class EdgeScoreBlender(EdgeScore):
	"""
	Blends two attribute vectors, the value is chosen depending on the supplied bool vector

	Parameters
	----------
	G : networkit.Graph
		The graph for which the attribute shall be blended
	attribute0 : vector[double]
		The first attribute (chosen for selection[eid] == false)
	attribute1 : vector[double]
		The second attribute (chosen for selection[eid] == true)
	selection : vector[bool]
		The selection vector
	"""
	cdef vector[double] _attribute0
	cdef vector[double] _attribute1
	cdef vector[bool_t] _selection

	def __cinit__(self, Graph G not None, vector[double] attribute0, vector[double] attribute1, vector[bool_t] selection):
		self._G = G
		self._attribute0 = move(attribute0)
		self._attribute1 = move(attribute1)
		self._selection = move(selection)

		self._this = new _EdgeScoreBlender(G._this, self._attribute0, self._attribute1, self._selection)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/edgescores/GeometricMeanScore.hpp>":

	cdef cppclass _GeometricMeanScore "NetworKit::GeometricMeanScore"(_EdgeScore[double]):
		_GeometricMeanScore(const _Graph& G, const vector[double]& a) except +

cdef class GeometricMeanScore(EdgeScore):
	"""
	Normalizes the given edge attribute by the geometric average of the sum of the attributes of the incident edges of the incident nodes.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	a : vector[double]
		Edge attribute that shall be normalized.
	"""
	cdef vector[double] _attribute

	def __cinit__(self, Graph G, vector[double] attribute):
		self._G = G
		self._attribute = attribute
		self._this = new _GeometricMeanScore(G._this, self._attribute)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/edgescores/EdgeScoreAsWeight.hpp>":

	cdef cppclass _EdgeScoreAsWeight "NetworKit::EdgeScoreAsWeight":
		_EdgeScoreAsWeight(const _Graph& G, const vector[double]& score, bool_t squared, edgeweight offset, edgeweight factor) except +
		_Graph calculate() except +

cdef class EdgeScoreAsWeight:
	"""
	Assigns an edge score as edge weight of a graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to assign edge weights to.
	score : vector[double]
		The input edge score.
	squared : bool
		Edge weights will be squared if set to True.
	offset : edgeweight
		This offset will be added to each edge weight.
	factor : edgeweight
		Each edge weight will be multiplied by this factor.
	"""

	cdef _EdgeScoreAsWeight* _this
	cdef Graph _G
	cdef vector[double] _score

	def __cinit__(self, Graph G, vector[double] score, bool_t squared, edgeweight offset, edgeweight factor):
		self._G = G
		self._score = score
		self._this = new _EdgeScoreAsWeight(G._this, self._score, squared, offset, factor)

	def __dealloc__(self):
		if self._this is not NULL:
			del self._this
			self._this = NULL

	def getWeightedGraph(self):
		"""
		Returns
		-------
		networkit.Graph
			The weighted result graph.
		"""
		return Graph(0).setThis(self._this.calculate())

# Module: sparsification

cdef extern from "<networkit/sparsification/SimmelianOverlapScore.hpp>":

	cdef cppclass _SimmelianOverlapScore "NetworKit::SimmelianOverlapScore"(_EdgeScore[double]):
		_SimmelianOverlapScore(const _Graph& G, const vector[count]& triangles, count maxRank) except +

cdef class SimmelianOverlapScore(EdgeScore):
	cdef vector[count] _triangles

	"""
	An implementation of the parametric variant of Simmelian Backbones. Calculates
	for each edge the minimum parameter value such that the edge is still contained in
	the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Simmelian Backbone algorithm to.
	triangles : vector[count]
		Previously calculated edge triangle counts on G.
	maxRank: count
		maximum rank that is considered for overlap calculation.
	"""
	def __cinit__(self, Graph G, vector[count] triangles, count maxRank):
		self._G = G
		self._triangles = triangles
		self._this = new _SimmelianOverlapScore(G._this, self._triangles, maxRank)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/edgescores/PrefixJaccardScore.hpp>":

	cdef cppclass _PrefixJaccardScore "NetworKit::PrefixJaccardScore<double>"(_EdgeScore[double]):
		_PrefixJaccardScore(const _Graph& G, const vector[double]& a) except +

cdef class PrefixJaccardScore(EdgeScore):
	cdef vector[double] _attribute

	def __cinit__(self, Graph G, vector[double] attribute):
		self._G = G
		self._attribute = attribute
		self._this = new _PrefixJaccardScore(G._this, self._attribute)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/MultiscaleScore.hpp>":

	cdef cppclass _MultiscaleScore "NetworKit::MultiscaleScore"(_EdgeScore[double]):
		_MultiscaleScore(const _Graph& G, const vector[double]& a) except +

cdef class MultiscaleScore(EdgeScore):
	"""
	An implementation of the Multiscale Backbone. Calculates for each edge the minimum
	parameter value such that the edge is still contained in the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Multiscale algorithm to.
	attribute : vector[double]
		The edge attribute the Multiscale algorithm is to be applied to.
	"""

	cdef vector[double] _attribute

	def __cinit__(self, Graph G, vector[double] attribute):
		self._G = G
		self._attribute = attribute
		self._this = new _MultiscaleScore(G._this, self._attribute)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/RandomEdgeScore.hpp>":

	cdef cppclass _RandomEdgeScore "NetworKit::RandomEdgeScore"(_EdgeScore[double]):
		_RandomEdgeScore(const _Graph& G) except +

cdef class RandomEdgeScore(EdgeScore):
	"""
	Generates a random edge attribute. Each edge is assigned a random value in [0,1].

	Parameters
	----------
	G : networkit.Graph
		The graph to calculate the Random Edge attribute for.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _RandomEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/LocalSimilarityScore.hpp>":

	cdef cppclass _LocalSimilarityScore "NetworKit::LocalSimilarityScore"(_EdgeScore[double]):
		_LocalSimilarityScore(const _Graph& G, const vector[count]& triangles) except +

cdef class LocalSimilarityScore(EdgeScore):
	"""
	An implementation of the Local Simlarity sparsification approach.
	This attributizer calculates for each edge the maximum parameter value
	such that the edge is still contained in the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Local Similarity algorithm to.
	triangles : vector[count]
		Previously calculated edge triangle counts.
	"""
	cdef vector[count] _triangles

	def __cinit__(self, Graph G, vector[count] triangles):
		self._G = G
		self._triangles = triangles
		self._this = new _LocalSimilarityScore(G._this, self._triangles)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/ForestFireScore.hpp>":

	cdef cppclass _ForestFireScore "NetworKit::ForestFireScore"(_EdgeScore[double]):
		_ForestFireScore(const _Graph& G, double pf, double tebr) except +

cdef class ForestFireScore(EdgeScore):
	"""
	A variant of the Forest Fire sparsification approach that is based on random walks.
	This attributizer calculates for each edge the minimum parameter value
	such that the edge is still contained in the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Forest Fire algorithm to.
	pf : double
		The probability for neighbor nodes to get burned aswell.
	tebr : double
		The Forest Fire will burn until tebr * numberOfEdges edges have been burnt.
	"""

	def __cinit__(self, Graph G, double pf, double tebr):
		self._G = G
		self._this = new _ForestFireScore(G._this, pf, tebr)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/LocalDegreeScore.hpp>":

	cdef cppclass _LocalDegreeScore "NetworKit::LocalDegreeScore"(_EdgeScore[double]):
		_LocalDegreeScore(const _Graph& G) except +

cdef class LocalDegreeScore(EdgeScore):
	"""
	The LocalDegree sparsification approach is based on the idea of hub nodes.
	This attributizer calculates for each edge the maximum parameter value
	such that the edge is still contained in the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Local Degree  algorithm to.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _LocalDegreeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/RandomNodeEdgeScore.hpp>":

	cdef cppclass _RandomNodeEdgeScore "NetworKit::RandomNodeEdgeScore"(_EdgeScore[double]):
		_RandomNodeEdgeScore(const _Graph& G) except +

cdef class RandomNodeEdgeScore(EdgeScore):
	"""
	Random Edge sampling. This attributizer returns edge attributes where
	each value is selected uniformly at random from [0,1].

	Parameters
	----------
	G : networkit.Graph
		The graph to calculate the Random Edge attribute for.
	"""
	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _RandomNodeEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return True

ctypedef fused DoubleInt:
	int
	double

cdef extern from "<networkit/sparsification/LocalFilterScore.hpp>":

	cdef cppclass _LocalFilterScoreDouble "NetworKit::LocalFilterScore<double>"(_EdgeScore[double]):
		_LocalFilterScoreDouble(const _Graph& G, const vector[double]& a, bool_t logarithmic) except +

	cdef cppclass _LocalFilterScoreInt "NetworKit::LocalFilterScore<int>"(_EdgeScore[count]):
		_LocalFilterScoreInt(const _Graph& G, const vector[double]& a, bool_t logarithmic) except +

cdef class LocalFilterScore(EdgeScore):
	"""
	Local filtering edge scoring. Edges with high score are more important.

	Edges are ranked locally, the top d^e (logarithmic, default) or 1+e*(d-1) edges (non-logarithmic) are kept.
	For equal attribute values, neighbors of low degree are preferred.
	If bothRequired is set (default: false), both neighbors need to indicate that they want to keep the edge.

	Parameters
	----------
	G : networkit.Graph
		The input graph
	a : list
		The input attribute according to which the edges shall be fitlered locally.
	logarithmic : bool
		If the score shall be logarithmic in the rank (then d^e edges are kept). Linear otherwise.
	"""
	cdef vector[double] _a

	def __cinit__(self, Graph G, vector[double] a, bool_t logarithmic = True):
		self._G = G
		self._a = a
		self._this = new _LocalFilterScoreDouble(G._this, self._a, logarithmic)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/ChanceCorrectedTriangleScore.hpp>":

	cdef cppclass _ChanceCorrectedTriangleScore "NetworKit::ChanceCorrectedTriangleScore"(_EdgeScore[double]):
		_ChanceCorrectedTriangleScore(const _Graph& G, const vector[count]& triangles) except +

cdef class ChanceCorrectedTriangleScore(EdgeScore):
	"""
	Divide the number of triangles per edge by the expected number of triangles given a random edge distribution.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	triangles : vector[count]
		Triangle count.
	"""
	cdef vector[count] _triangles

	def __cinit__(self, Graph G, vector[count] triangles):
		self._G = G
		self._triangles = triangles
		self._this = new _ChanceCorrectedTriangleScore(G._this, self._triangles)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/SCANStructuralSimilarityScore.hpp>":

	cdef cppclass _SCANStructuralSimilarityScore "NetworKit::SCANStructuralSimilarityScore"(_EdgeScore[double]):
		_SCANStructuralSimilarityScore(_Graph G, const vector[count]& triangles) except +

cdef class SCANStructuralSimilarityScore(EdgeScore):
	"""
	An implementation of the SCANStructuralSimilarityScore algorithm.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Local Similarity algorithm to.
	triangles : vector[count]
		Previously calculated edge triangle counts.
	"""
	cdef vector[count] _triangles

	def __cinit__(self, Graph G, vector[count] triangles):
		self._G = G
		self._triangles = triangles
		self._this = new _SCANStructuralSimilarityScore(G._this, self._triangles)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/GlobalThresholdFilter.hpp>":

	cdef cppclass _GlobalThresholdFilter "NetworKit::GlobalThresholdFilter":
		_GlobalThresholdFilter(const _Graph& G, const vector[double]& a, double alpha, bool_t above) except +
		_Graph calculate() except +

cdef class GlobalThresholdFilter:
	"""
	Calculates a sparsified graph by filtering globally using a constant threshold value
	and a given edge attribute.

	Parameters
	----------
	G : networkit.Graph
		The graph to sparsify.
	attribute : vector[double]
		The edge attribute to consider for filtering.
	e : double
		Threshold value.
	above : bool
		If set to True (False), all edges with an attribute value equal to or above (below)
		will be kept in the sparsified graph.
	"""
	cdef _GlobalThresholdFilter* _this
	cdef Graph _G
	cdef vector[double] _attribute

	def __cinit__(self, Graph G not None, vector[double] attribute, double e, bool_t above):
		self._G = G
		self._attribute = attribute
		self._this = new _GlobalThresholdFilter(G._this, self._attribute, e, above)

	def __dealloc__(self):
		del self._this

	def calculate(self):
		return Graph().setThis(self._this.calculate())


cdef extern from "<networkit/matching/Matcher.hpp>":

	cdef cppclass _Matcher "NetworKit::Matcher"(_Algorithm):
		_Matcher(const _Graph _G) except +
		_Matching getMatching() except +

cdef class Matcher(Algorithm):
	""" Abstract base class for matching algorithms """
	cdef Graph G

	def __init__(self, *args, **namedargs):
		if type(self) == Matcher:
			raise RuntimeError("Instantiation of abstract base class")

	def getMatching(self):
		"""  Returns the matching.

		Returns
		-------
		Matching
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return Matching().setThis((<_Matcher*>(self._this)).getMatching())


cdef extern from "<networkit/matching/PathGrowingMatcher.hpp>":

	cdef cppclass _PathGrowingMatcher "NetworKit::PathGrowingMatcher"(_Matcher):
		_PathGrowingMatcher(_Graph) except +
		_PathGrowingMatcher(_Graph, vector[double]) except +

cdef class PathGrowingMatcher(Matcher):
	"""
	Path growing matching algorithm as described by  Hougardy and Drake.
	Computes an approximate maximum weight matching with guarantee 1/2.
	"""
	def __cinit__(self, Graph G not None, edgeScores=None):
		self.G = G
		if edgeScores:
			self._this = new _PathGrowingMatcher(G._this, edgeScores)
		else:
			self._this = new _PathGrowingMatcher(G._this)

# profiling

def ranked(sample):
	"""
		Given a list of numbers, this function computes the rank of each value
		and returns a list of ranks where result[i] is the rank of
		the i-th element in the given sample.
		Currently used in profiling.stat.
	"""
	cdef vector[pair[double, count]] helper = vector[pair[double, count]](len(sample))
	cdef vector[double] result = vector[double](len(sample), 0)
	for i in range(len(sample)):
		helper[i] = <pair[double, count]?>(sample[i], i)
	sort(helper.begin(), helper.end())
	cdef double value = helper[0].first
	cdef double summ = 0.
	cdef count length = 0
	for i in range(len(sample)):
		if value == helper[i].first:
			summ += (i+1)
			length += 1
		else:
			summ /= length
			for j in range(length):
				result[helper[i-j-1].second] = summ
			value = helper[i].first
			summ = i+1
			length = 1
	summ /= length
	for j in range(length):
		result[helper[len(sample)-j-1].second] = summ
	return result

def sort2(sample):
	"""
		Sorts a given list of numbers.
		Currently used as profiling.stat.sorted.
	"""
	cdef vector[double] result = <vector[double]?>sample
	sort(result.begin(),result.end())
	return result

# stats

def gini(values):
	"""
	Computes the Gini coefficient for the distribution given as a list of values.
	"""
	sorted_list = sorted(values)
	height, area = 0, 0
	for value in sorted_list:
		height += value
		area += height - value / 2.
	fair_area = height * len(values) / 2
	return (fair_area - area) / fair_area


# simulation
cdef extern from "<networkit/simulation/EpidemicSimulationSEIR.hpp>":

	cdef cppclass _EpidemicSimulationSEIR "NetworKit::EpidemicSimulationSEIR" (_Algorithm):
		_EpidemicSimulationSEIR(_Graph, count, double, count, count, node) except +
		vector[vector[count]] getData() except +

cdef class EpidemicSimulationSEIR(Algorithm):
	"""
 	Parameters
 	----------
 	G : networkit.Graph
 		The graph.
 	tMax : count
 		max. number of timesteps
	transP : double
		transmission probability
	eTime : count
		exposed time
	iTime : count
		infectious time
	zero : node
		starting node
	"""
	cdef Graph G
	def __cinit__(self, Graph G, count tMax, double transP=0.5, count eTime=2, count iTime=7, node zero=none):
		self.G = G
		self._this = new _EpidemicSimulationSEIR(G._this, tMax, transP, eTime, iTime, zero)
	def getData(self):
		return pandas.DataFrame((<_EpidemicSimulationSEIR*>(self._this)).getData(), columns=["zero", "time", "state", "count"])



## Module: viz


cdef extern from "<networkit/viz/GraphLayoutAlgorithm.hpp>":

	cdef cppclass _GraphLayoutAlgorithm "NetworKit::GraphLayoutAlgorithm"[T]:
		_GraphLayoutAlgorithm(_Graph, count) except +
		count numEdgeCrossings() except +
		vector[Point[double]] getCoordinates() except +
		bool_t writeGraphToGML(string path) except +
		bool_t writeKinemage(string path) except +

cdef class GraphLayoutAlgorithm:

	"""Abstract base class for graph drawing algorithms"""

	cdef _GraphLayoutAlgorithm[double] *_this
	cdef Graph _G

	def __init__(self, *args, **kwargs):
		if type(self) == GraphLayoutAlgorithm:
			raise RuntimeError("Error, you may not use GraphLayoutAlgorithm directly, use a sub-class instead")

	def __dealloc__(self):
		self._G = None # just to be sure the graph is deleted

	def numEdgeCrossings(self):
		""" Computes approximation (in parallel) of the Spanning Edge Centrality. """
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.numEdgeCrossings()

	def getCoordinates(self):
		""" Computes approximation (in parallel) of the Spanning Edge Centrality. """
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		cdef pair[double, double] pr = pair[double, double](0, 0)
		pointCoord = self._this.getCoordinates()
		cdef vector[pair[double, double]] pairCoord = vector[pair[double, double]]()
		for pt in pointCoord:
			pr = pair[double, double](pt[0], pt[1])
			pairCoord.push_back(pr)
		return pairCoord

	def writeGraphToGML(self, path):
		"""Writes the graph and its layout to a .gml file at the specified path
	path: string
		Path where the graph file should be created"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.writeGraphToGML(stdstring(path))

	def writeKinemage(self, string path):
		"""Writes the graph and its layout to a file at the specified path
			path: string
		Path where the graph file should be created"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.writeKinemage(stdstring(path))



cdef extern from "<networkit/viz/MaxentStress.hpp>" namespace "NetworKit":

	enum _GraphDistance "NetworKit::MaxentStress::GraphDistance":
		EDGE_WEIGHT,
		ALGEBRAIC_DISTANCE

cdef extern from "<networkit/viz/MaxentStress.hpp>" namespace "NetworKit":

	enum _LinearSolverType "NetworKit::MaxentStress::LinearSolverType":
		LAMG,
		CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER,
		CONJUGATE_GRADIENT_DIAGONAL_PRECONDITIONER


cdef extern from "<networkit/viz/MaxentStress.hpp>":

	cdef cppclass _MaxentStress "NetworKit::MaxentStress" (_GraphLayoutAlgorithm[double]):
		_MaxentStress(_Graph G, count dim, count k, double tolerance, _LinearSolverType linearSolverType, bool_t fastComputation, _GraphDistance graphDistance) except +
		_MaxentStress(_Graph G, count dim, const vector[Point[double]] coordinates, count k, double tolerance, _LinearSolverType linearSolverType, bool_t fastComputation, _GraphDistance graphDistance) except +
		void run() except +
		void scaleLayout() except +
		double computeScalingFactor() except +
		double fullStressMeasure() except +
		double maxentMeasure() except +
		double meanDistanceError() except +
		double ldme() except +
		void setQ(double q) except +
		void setAlpha(double alpha) except +
		void setAlphaReduction(double alphaReduction) except +
		void setFinalAlpha(double finalAlpha) except +
		void setConvergenceThreshold(double convThreshold) except +
		double getRhs() except +
		double getApproxEntropyTerm() except +
		double getSolveTime() except +


cdef class MaxentStress (GraphLayoutAlgorithm):

	"""
	Implementation of MaxentStress by Gansner et al. using a Laplacian system solver.
  	@see Gansner, Emden R., Yifan Hu, and Steve North. "A maxent-stress model for graph layout."
	Visualization and Computer Graphics, IEEE Transactions on 19, no. 6 (2013): 927-940.

	Parameters
	----------
	G : networkit.Graph
		The graph to be handled. Should be connected, otherwise the run() and runAlgo() methods will fail.
	dim: count
		Number of dimensions.
	count: k
	coordinates: vector[pair[double, double]]
		The coordinates we want to draw in.
	tolerance: double
		The tolerance we want our solver to have.
	linearSolverType: _LinearSolverType
		The type of linear solver we wish to use.
	fastComputation: bool
		Decides whether or not slightly faster computation should be employed, leading to slightly worse results.
	graphDistance: _GraphDistance
		Decides what type of graph distance should be utilised.
	"""

	LAMG = 0
	CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER = 1
	CONJUGATE_GRADIENT_DIAGONAL_PRECONDITIONER = 2
	EDGE_WEIGHT = 0
	ALGEBRAIC_DISTANCE = 1

	def __cinit__(self, Graph G, count dim, count k, vector[pair[double, double]] coordinates = [], double tolerance = 1e-5, _LinearSolverType linearSolverType = LAMG, bool_t fastComputation = False, _GraphDistance graphDistance = EDGE_WEIGHT):
		cdef Point[double] p = Point[double](0, 0)
		cdef vector[Point[double]] pointCoordinates = vector[Point[double]]()

		for pr in coordinates:
			p = Point[double](pr.first, pr.second)
			pointCoordinates.push_back(p)

		if (coordinates.size() != 0):
			self._this = new _MaxentStress(G._this, dim, pointCoordinates, k, tolerance, linearSolverType, fastComputation, graphDistance)
		else:
			self._this = new _MaxentStress(G._this, dim, k, tolerance, linearSolverType, fastComputation, graphDistance)

	def __dealloc__(self):
		del self._this

	def run(self):
		"""Approximates a graph layout with the maxent-stress algorithm"""
		(<_MaxentStress*>(self._this)).run()
		return self

	def scaleLayout(self):
		"""Scale the layout computed by run() by a scalar s to minimize \sum_{u,v \in V} w_{uv} (s ||x_u - x_v|| - d_{uv}||)^2"""
		(<_MaxentStress*>(self._this)).scaleLayout()
		return self

	def computeScalingFactor(self):
		"""Computes a scalar s s.t. \sum_{u,v \in V} w_{uv} (s ||x_u - x_v|| - d_{uv}||)^2 is minimized"""
		return (<_MaxentStress*>(self._this)).computeScalingFactor()

	def fullStressMeasure(self):
		"""Computes the full stress measure of the computed layout with run()"""
		return (<_MaxentStress*>(self._this)).fullStressMeasure()

	def maxentMeasure(self):
		"""Computes the maxent stress measure for the computed layout with run()"""
		return (<_MaxentStress*>(self._this)).maxentMeasure()

	def meanDistanceError(self):
		"""Computes mean distance error"""
		return (<_MaxentStress*>(self._this)).meanDistanceError()

	def ldme(self):
		"""Computes the ldme"""
		return (<_MaxentStress*>(self._this)).ldme()

	def setQ(self, double q):
		(<_MaxentStress*>(self._this)).setQ(q)
		return self

	def setAlpha(self, double alpha):
		(<_MaxentStress*>(self._this)).setAlpha(alpha)
		return self

	def setAlphaReduction(self, double alphaReduction):
		(<_MaxentStress*>(self._this)).setAlphaReduction(alphaReduction)
		return self

	def setFinalAlpha(self, double finalAlpha):
		(<_MaxentStress*>(self._this)).setFinalAlpha(finalAlpha)
		return self

	def setConvergenceThreshold(self, double convThreshold):
		(<_MaxentStress*>(self._this)).setConvergenceThreshold(convThreshold)
		return self

	def getRhs(self):
		return (<_MaxentStress*>(self._this)).getRhs()

	def getApproxEntropyTerm(self):
		return (<_MaxentStress*>(self._this)).getApproxEntropyTerm()

	def getSolveTime(self):
		return (<_MaxentStress*>(self._this)).getSolveTime()





cdef extern from "<networkit/viz/PivotMDS.hpp>":

	cdef cppclass _PivotMDS "NetworKit::PivotMDS" (_GraphLayoutAlgorithm[double]):
				_PivotMDS(_Graph G, count dim, count numberOfPivots) except +
				void run() except +


cdef class PivotMDS (GraphLayoutAlgorithm):

	"""
	Implementation of PivotMDS proposed by Brandes and Pich.

	Parameters
	----------

	G: networkit.Graph
		The graph to be handled by the algorithm.

	dim: count
		Number of dimensions.

	numberOfPivots: count
		Number of pivots for the algorithm.

	"""

	def __cinit__(self, Graph G, count dim, count numberOfPivots):
		self._this = new _PivotMDS(G._this, dim, numberOfPivots)

	def __dealloc__(self):
		del self._this

	def run(self):
		"""Constructs a PivotMDS object for the given @a graph. The algorithm should embed the graph in @a dim dimensions using @a numberOfPivots pivots."""
		(<_PivotMDS*>(self._this)).run()
		return self


# Module: randomization

cdef extern from "<networkit/randomization/GlobalCurveball.hpp>":

	cdef cppclass _GlobalCurveball "NetworKit::GlobalCurveball"(_Algorithm):
		_GlobalCurveball(_Graph, count, bool_t, bool_t) except +
		_Graph getGraph() except +

cdef class GlobalCurveball(Algorithm):
	"""
	Implementation of EM-GCB proposed in "Parallel and I/O-efficient
	Randomisation of Massive Networks using Global Curveball Trades",
	Carstens et al., ESA 2018.

	The algorithm perturbs an unweighted input graph, by iteratively
	randomizing the neighbourhoods of node pairs. For a large number
	of global trades this process is shown to produce an uniform sample
	from the set of all graphs with the same degree sequence as the input
	graph.

	If you do not want to explicitly control the trade sequence,
	we recommend using GlobalCurveball rather than Curveball since
	GlobalCurveball is typically faster and exhibits a smaller memory
	footprint.

	Parameters
	----------

	G : networkit.Graph
		The graph to be randomized. For a given degree sequence, e.g.
		generators.HavelHakimi can be used to obtain this graph.

	number_of_global_rounds:
		Number of global rounds to carry out. The runtime scales
		asymptotically linearly in this parameter. Default: 20,
		which yields good results experimentally (see Paper).

	allowSelfLoops:
		Has to be False for undirected graphs. For directed graphs
		the randomization Markov chain is only irreducible if self loops
		are allows. If they are forbidden, the degreePreservingShuffle
		proprocessing has to be enabled. Otherwhise, not all topologies
		can be produced.

	degreePreservingShufflePreprocessing:
		Execute the DegreePreservingShuffle algorithm before executing
		Global Curveball. It's more efficient than manually invoking
		the algorithm.

	Warning
	-------
	For directed graphs at least one of allowSelfLoops or
	degreePreservingShufflePreprocessing should be set; for more details
	refer to "Switching edges to randomize networks: what goes wrong
	and how to fix it" by C. J. Carstens K. J. Horadam
	"""
	def __cinit__(self, G, number_of_global_rounds = 20, allowSelfLoops = False, degreePreservingShufflePreprocessing = True):
		if isinstance(G, Graph):
			self._this = new _GlobalCurveball((<Graph>G)._this, number_of_global_rounds, allowSelfLoops, degreePreservingShufflePreprocessing)
		else:
			raise RuntimeError("Parameter G has to be a graph")

	"""

	Get randomized graph after invocation of run().

	"""

	def getGraph(self):
		return Graph().setThis((<_GlobalCurveball*>self._this).getGraph())


cdef extern from "<networkit/randomization/CurveballUniformTradeGenerator.hpp>":

	cdef cppclass _CurveballUniformTradeGenerator "NetworKit::CurveballUniformTradeGenerator":
		_CurveballUniformTradeGenerator(count runLength, count numNodes) except +
		vector[pair[node, node]] generate() nogil except +

cdef class CurveballUniformTradeGenerator:

	"""
	Generates a trade sequence consisting of num_trades many single trades.
	Each trade contains two different node indices drawn uniformly at random
	from the interval [0, num_nodes).

	Parameters
	----------

	num_trades:
	   Number of trades to generate.

	num_nodes:
	   Number of node indices to draw from

	"""
	cdef _CurveballUniformTradeGenerator *_this

	def __cinit__(self, count num_trades, count num_nodes):
		self._this = new _CurveballUniformTradeGenerator(num_trades, num_nodes)

	def __dealloc__(self):
		del self._this

	def generate(self):
		return self._this.generate()

cdef extern from "<networkit/randomization/CurveballGlobalTradeGenerator.hpp>":

	cdef cppclass _CurveballGlobalTradeGenerator "NetworKit::CurveballGlobalTradeGenerator":
		_CurveballGlobalTradeGenerator(count runLength, count numNodes) except +
		vector[pair[node, node]] generate() nogil except +

cdef class CurveballGlobalTradeGenerator:

	"""
	Generates a trade sequence consisting of num_global_trades global trades
	targeting node ids from the range [0, num_nods).

	If you are only using this generator, consider using the GlobalCurveball
	algorithm directly as it has a better performance / memory footprint.

	Parameters
	----------

	num_global_trades:
	   Number of global trades to generate (i.e. the resulting sequence contains
	   num_global_trades * floor(num_nodes / 2) trades)

	num_nodes:
	   Number of node indices to draw from

	"""
	cdef _CurveballGlobalTradeGenerator *_this

	def __cinit__(self, count num_global_trades, count num_nodes):
		self._this = new _CurveballGlobalTradeGenerator(num_global_trades, num_nodes)

	def __dealloc__(self):
		del self._this

	def generate(self):
		return self._this.generate()

cdef extern from "<networkit/randomization/Curveball.hpp>":

	cdef cppclass _Curveball "NetworKit::Curveball"(_Algorithm):
		_Curveball(_Graph) except +
		void run(vector[pair[node, node]] trades) nogil except +
		_Graph getGraph() except +
		vector[pair[node, node]] getEdges() except +
		count getNumberOfAffectedEdges() except +

cdef class Curveball(Algorithm):
	"""
	Implementation of IM-CB proposed in "Parallel and I/O-efficient
	Randomisation of Massive Networks using Global Curveball Trades",
	Carstens et al., ESA 2018.

	The algorithm perturbs an undirected and unweighted input graph,
	by iteratively randomizing the neighbourhoods of node pairs. For
	a large number of trades this process is shown to produce an
	uniform sample from the set of all graphs with the same degree
	sequence as the input graph.

	If you do not want to explicitly control the trade sequence,
	we recommend using GlobalCurveball rather than Curveball since
	GlobalCurveball is typically faster and exhibits a smaller memory
	footprint.

   Observe that this algorithm does not support the run() method,
   since it requires the trade sequence to be passed. It is possible
   to invoke run(trades) several times, e.g. to reduce the memory
   footprint which increases linearly with the number of trades
   performed in a run.

	Parameters
	----------

	G : networkit.Graph
		The graph to be randomized. For a given degree sequence, e.g.
		generators.HavelHakimi can be used to obtain this graph.

	"""
	def __cinit__(self, G):
		if isinstance(G, Graph):
			self._this = new _Curveball((<Graph>G)._this)
		else:
			raise RuntimeError("Parameter G has to be a graph")

	def run(self, vector[pair[node, node]] trades):
		with nogil:
			(<_Curveball*>(self._this)).run(trades)
		return self

	def getGraph(self):
		return Graph().setThis((<_Curveball*>self._this).getGraph())

	def getNumberOfAffectedEdges(self):
		return (<_Curveball*>(self._this)).getNumberOfAffectedEdges()


cdef extern from "../include/networkit/randomization/DegreePreservingShuffle.hpp":
	cdef cppclass _DegreePreservingShuffle "NetworKit::DegreePreservingShuffle"(_Algorithm):
		_DegreePreservingShuffle(_Graph) except +
		_Graph getGraph() except +
		vector[node] getPermutation() except +

cdef class DegreePreservingShuffle(Algorithm):
	"""
	Implementation of the preprocessing step proposed in
	"Smaller Universes for Uniform Sampling of 0,1-matrices with fixed row and column sums"
	by Annabell Berger, Corrie Jacobien Carstens [https://arxiv.org/abs/1803.02624]

	The algorithms randomizes a graph without changing its topology simply
	by renaming nodes. For any degree d (in case of an directed graph it's a degree pair)
	consider the set X_d of node ids which have this degree. Then shuffle the ids in X_d.

	Hence the algorithm satisfies: For all x in Ginput:
	 i)  Ginput.degreeIn(x) = Goutput.degreeIn(x)
	 ii) Ginput.degreeOut(x) = Goutput.degreeOut(x)

	The authors argue that applying this preprocessing step before executing (Global)Curveball
	leads to a faster mixing time. If you want to use it as a preprocessing step to GlobalCurveball,
	it's more efficient to set degreePreservingShufflePreprocessing in GlobalCurveball's constructor.

	Parameters
	----------

	G : networkit.Graph
		The graph to be randomized. For a given degree sequence, e.g.
		generators.HavelHakimi can be used to obtain this graph.

	"""
	def __cinit__(self, G):
		if isinstance(G, Graph):
			self._this = new _DegreePreservingShuffle((<Graph>G)._this)
		else:
			raise RuntimeError("Parameter G has to be a graph")

	def getGraph(self):
		return Graph().setThis((<_DegreePreservingShuffle*>self._this).getGraph())

	def getPermutation(self):
		return (<_DegreePreservingShuffle*>(self._this)).getPermutation()
