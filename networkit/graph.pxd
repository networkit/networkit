# distutils: language=c++

from cython.operator import dereference, preincrement

from libc.stdint cimport uint64_t

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.unordered_set cimport unordered_set

ctypedef uint64_t edgeid
ctypedef uint64_t index
ctypedef uint64_t count
ctypedef index node
ctypedef double edgeweight

from .base cimport _Algorithm
from .base cimport Algorithm

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Graph move( _Graph t ) nogil
	vector[double] move(vector[double])

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

cdef extern from "<networkit/graph/Graph.hpp>":

	cdef struct Edge "NetworKit::Edge":
		node u
		node v

	cdef struct WeightedEdge "NetworKit::WeightedEdge":
		node u
		node v
		edgeweight weight

	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(count, bool_t, bool_t) except +
		_Graph(const _Graph& other) except +
		_Graph(const _Graph& other, bool_t weighted, bool_t directed) except +
		void indexEdges(bool_t) except +
		bool_t hasEdgeIds() except +
		edgeid edgeId(node, node) except +
		count numberOfNodes() except +
		count numberOfEdges() except +
		index upperNodeIdBound() except +
		index upperEdgeIdBound() except +
		count degree(node u) except +
		count degreeIn(node u) except +
		count degreeOut(node u) except +
		double weightedDegree(node u, bool_t) except +
		double weightedDegreeIn(node u, bool_t) except +
		bool_t isIsolated(node u) except +
		node addNode() except +
		node addNodes(node) except +
		void removeNode(node u) except +
		bool_t hasNode(node u) except +
		void restoreNode(node u) except +
		void addEdge(node u, node v, edgeweight w) except +
		void setWeight(node u, node v, edgeweight w) except +
		void increaseWeight(node u, node v, edgeweight w) except +
		void removeEdge(node u, node v) except +
		void removeAllEdges() except +
		void removeSelfLoops() except +
		void removeMultiEdges() except +
		void swapEdge(node s1, node t1, node s2, node t2) except +
		void compactEdges() except +
		void sortEdges() except +
		bool_t hasEdge(node u, node v) except +
		edgeweight weight(node u, node v) except +
		void forEdges[Callback](Callback c) except +
		void forNodes[Callback](Callback c) except +
		void forNodePairs[Callback](Callback c) except +
		void forNodesInRandomOrder[Callback](Callback c) except +
		void forEdgesOf[Callback](node u, Callback c) except +
		void forInEdgesOf[Callback](node u, Callback c) except +
		bool_t isWeighted() except +
		bool_t isDirected() except +
		edgeweight totalEdgeWeight() except +
		count numberOfSelfLoops() except +
		bool_t checkConsistency() except +
		_NodeRange nodeRange() except +
		_EdgeRange edgeRange() except +
		_EdgeWeightRange edgeWeightRange() except +
		_OutNeighborRange neighborRange(node u) except +
		_InNeighborRange inNeighborRange(node u) except +

cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _NodeIterator "NetworKit::Graph::NodeIterator":
		_NodeIterator operator++() except +
		_NodeIterator operator++(int) except +
		bool_t operator!=(const _NodeIterator) except +
		node operator*() except +

cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _NodeRange "NetworKit::Graph::NodeRange":
		_NodeIterator begin() except +
		_NodeIterator end() except +

cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _EdgeWeightIterator "NetworKit::Graph::EdgeWeightIterator":
		_EdgeWeightIterator operator++() except +
		_EdgeWeightIterator operator++(int) except +
		bool_t operator!=(const _EdgeWeightIterator) except +
		WeightedEdge operator*() except +

cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _EdgeWeightRange "NetworKit::Graph::EdgeWeightRange":
		_EdgeWeightIterator begin() except +
		_EdgeWeightIterator end() except +

cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _EdgeIterator "NetworKit::Graph::EdgeIterator":
		_EdgeIterator operator++() except +
		_EdgeIterator operator++(int) except +
		bool_t operator!=(const _EdgeIterator) except +
		Edge operator*() except +

cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _EdgeRange "NetworKit::Graph::EdgeRange":
		_EdgeIterator begin() except +
		_EdgeIterator end() except +

cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _NeighborIterator "NetworKit::Graph::NeighborIterator":
		_NeighborIterator operator++() except +
		_NeighborIterator operator++(int) except +
		bool_t operator!=(const _NeighborIterator) except +
		node operator*() except +

cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _OutNeighborRange "NetworKit::Graph::OutNeighborRange":
		_NeighborIterator begin() except +
		_NeighborIterator end() except +

cdef extern from "<networkit/graph/Graph.hpp>":
	cdef cppclass _InNeighborRange "NetworKit::Graph::InNeighborRange":
		_NeighborIterator begin() except +
		_NeighborIterator end() except +

cdef class Graph:
	cdef _Graph _this
	cdef setThis(self, _Graph& other)

cdef extern from "<networkit/graph/SpanningForest.hpp>":

	cdef cppclass _SpanningForest "NetworKit::SpanningForest":
		_SpanningForest(_Graph) except +
		void run() nogil except +
		_Graph getForest() except +


cdef extern from "<networkit/graph/RandomMaximumSpanningForest.hpp>":

	cdef cppclass _RandomMaximumSpanningForest "NetworKit::RandomMaximumSpanningForest"(_Algorithm):
		_RandomMaximumSpanningForest(_Graph) except +
		_RandomMaximumSpanningForest(_Graph, vector[double]) except +
		_Graph getMSF(bool_t move) except +
		vector[bool_t] getAttribute(bool_t move) except +
		bool_t inMSF(edgeid eid) except +
		bool_t inMSF(node u, node v) except +

cdef class RandomMaximumSpanningForest(Algorithm):
	cdef vector[double] _attribute
	cdef Graph _G


cdef extern from "<networkit/graph/UnionMaximumSpanningForest.hpp>":

	cdef cppclass _UnionMaximumSpanningForest "NetworKit::UnionMaximumSpanningForest"(_Algorithm):
		_UnionMaximumSpanningForest(_Graph) except +
		_UnionMaximumSpanningForest(_Graph, vector[double]) except +
		_Graph getUMSF(bool_t move) except +
		vector[bool_t] getAttribute(bool_t move) except +
		bool_t inUMSF(edgeid eid) except +
		bool_t inUMSF(node u, node v) except +

cdef class UnionMaximumSpanningForest(Algorithm):
	cdef Graph _G

