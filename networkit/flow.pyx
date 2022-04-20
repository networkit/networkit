# distutils: language=c++

from libcpp.vector cimport vector

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport index, edgeid, node, edgeweight

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

none = _none

cdef extern from "<networkit/flow/EdmondsKarp.hpp>":

	cdef cppclass _EdmondsKarp "NetworKit::EdmondsKarp"(_Algorithm):
		_EdmondsKarp(const _Graph &graph, node source, node sink) except +
		edgeweight getMaxFlow() const
		vector[node] getSourceSet() except +
		edgeweight getFlow(node u, node v) except +
		edgeweight getFlow(edgeid eid) const
		vector[edgeweight] getFlowVector() except +

cdef class EdmondsKarp(Algorithm):
	"""
	EdmondsKarp(graph, source, sink)

	The EdmondsKarp class implements the maximum flow algorithm by Edmonds and Karp.

	Parameters
	----------
	graph : networkit.Graph
		The graph
	source : int
		The source node for the flow calculation
	sink : int
		The sink node for the flow calculation
	"""
	cdef Graph _graph

	def __cinit__(self, Graph graph not None, node source, node sink):
		self._graph = graph # store reference of graph for memory management, so the graph is not deallocated before this object
		self._this = new _EdmondsKarp(graph._this, source, sink)

	def getMaxFlow(self):
		"""
		getMaxFlow()

		Returns the value of the maximum flow from source to sink.

		Returns
		-------
		float
			The maximum flow value
		"""
		return (<_EdmondsKarp*>(self._this)).getMaxFlow()

	def getSourceSet(self):
		"""
		getSourceSet()

		Returns the set of the nodes on the source side of the flow/minimum cut.

		Returns
		-------
		list(int)
			The set of nodes that form the (smallest) source side of the flow/minimum cut.
		"""
		return (<_EdmondsKarp*>(self._this)).getSourceSet()

	def getFlow(self, node u, node v = none):
		"""
		getFlow(u, v = None)

		Get the flow value between two nodes u and v or an edge identified by the edge id u.
		Warning: The variant with two edge ids is linear in the degree of u.

		Parameters
		----------
		u : int
			The first node incident to the edge or the edge id.
		v : int, optional
			The second node incident to the edge (optional if edge id is specified). Default: None

		Returns
		-------
		float
			The flow on the specified edge.
		"""
		if v == none: # Assume that node and edge ids are the same type
			return (<_EdmondsKarp*>(self._this)).getFlow(u)
		else:
			return (<_EdmondsKarp*>(self._this)).getFlow(u, v)

	def getFlowVector(self):
		"""
		getFlowVector()

		Return a copy of the flow values of all edges.

		Returns
		-------
		list(float)
			The flow values of all edges indexed by edge id.
		"""
		return (<_EdmondsKarp*>(self._this)).getFlowVector()
