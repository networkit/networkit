# cython: language_level=3

# type imports
from libc.stdint cimport uint64_t

# the C++ standard library
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

# NetworKit typedefs
ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node
ctypedef double coordinate

from .base cimport _Algorithm
from .base cimport Algorithm
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition
from .dynamics cimport _GraphEvent, GraphEvent

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

