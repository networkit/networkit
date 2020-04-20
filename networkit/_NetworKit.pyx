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
