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

# === Dinic binding ===========================================================

cdef extern from "<networkit/flow/Dinic.hpp>":

	cdef cppclass _Dinic "NetworKit::Dinic"(_Algorithm):
		_Dinic(const _Graph &graph, node src, node dst) except +
		edgeweight getMaxFlow() except +


cdef class Dinic(Algorithm):
	"""
	Dinic(graph, source, sink)

	Computes maximum flow in a directed, weighted graph using Dinic's
	blocking-flow algorithm.

	Parameters
	----------
	graph : networkit.Graph
		Directed, weighted input graph.
	source : int
		Source node identifier.
	sink : int
		Target node identifier.
	"""
	cdef Graph _graph

	def __cinit__(self, Graph graph not None, node source, node sink):
		self._graph = graph
		self._this = new _Dinic(graph._this, source, sink)

	def getMaxFlow(self):
		"""
		getMaxFlow()

		Returns the computed maximum flow from source to sink.

		Returns
		-------
		float
			The maximum flow value.
		"""
		return (<_Dinic*>(self._this)).getMaxFlow()

# === SuccessiveShortestPath (Min-Cost Flow) binding ==========================

cdef extern from "<string_view>" namespace "std":
	cdef cppclass string_view:
		string_view() except +
		string_view(const char*, size_t) except +

# Bind EdgeDoubleAttribute via its fully qualified C++ type name.
# Note: this is a C++ 'using' alias (not a separately declared class),
# but it still names a concrete C++ type.
cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _EdgeDoubleAttribute "NetworKit::Graph::EdgeDoubleAttribute":
		double get(index i) const
		index size() const

cdef class EdgeDoubleAttribute:
	"""
	Thin Python wrapper for NetworKit::Graph::EdgeDoubleAttribute
	"""
	cdef _EdgeDoubleAttribute _this
	cdef object _owner  # keep the producing object alive (SSP instance)

	@staticmethod
	cdef EdgeDoubleAttribute _from_cpp(_EdgeDoubleAttribute a, object owner):
		cdef EdgeDoubleAttribute obj = EdgeDoubleAttribute.__new__(EdgeDoubleAttribute)
		obj._this = a
		obj._owner = owner
		return obj

	def size(self):
		return self._this.size()

	def get(self, index i):
		return self._this.get(i)

	def __getitem__(self, index i):
		return self._this.get(i)


cdef extern from "<networkit/flow/SuccessiveShortestPath.hpp>":

	cdef cppclass _SuccessiveShortestPathMinCostFlow "NetworKit::SuccessiveShortestPathMinCostFlow"(_Algorithm):
		_SuccessiveShortestPathMinCostFlow(const _Graph &G,
										   string_view capacityName,
										   string_view supplyName) except +
		double getTotalCost() const
		_EdgeDoubleAttribute getFlow() const


cdef class SuccessiveShortestPathMinCostFlow(Algorithm):
	"""
	SuccessiveShortestPathMinCostFlow(graph, capacityAttributeName, supplyAttributeName)
	"""
	cdef Graph _graph

	def __cinit__(self, Graph graph not None, capacityAttributeName, supplyAttributeName):
		self._graph = graph  # keep input graph alive (consistent with other bindings)

		cdef bytes cap_b = (<str>capacityAttributeName).encode("utf-8")
		cdef bytes sup_b = (<str>supplyAttributeName).encode("utf-8")

		cdef const char* cap_c = cap_b
		cdef const char* sup_c = sup_b

		cdef string_view cap_sv = string_view(cap_c, <size_t>len(cap_b))
		cdef string_view sup_sv = string_view(sup_c, <size_t>len(sup_b))

		self._this = new _SuccessiveShortestPathMinCostFlow(graph._this, cap_sv, sup_sv)

	def getTotalCost(self):
		return (<_SuccessiveShortestPathMinCostFlow*>(self._this)).getTotalCost()

	def getFlow(self):
		# Return the C++ EdgeDoubleAttribute (wrapped), exactly like the C++ API.
		cdef _EdgeDoubleAttribute a = (<_SuccessiveShortestPathMinCostFlow*>(self._this)).getFlow()
		return EdgeDoubleAttribute._from_cpp(a, self)
