# cython: language_level=3

#includes
# needed for collections.Iterable
from networkit.exceptions import ReducedFunctionalityWarning
import collections
import math
import os
import tempfile
import warnings

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

