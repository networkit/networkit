# distutils: language=c++

from libc.stdint cimport uint64_t
from cython.operator import dereference, preincrement

from libcpp cimport bool as bool_t
from libcpp.string cimport string
from libcpp.vector cimport vector

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node
ctypedef uint64_t edgeid
ctypedef double edgeweight

from .graph cimport _Graph, Graph
from .helpers import stdstring

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef cppclass TraversalEdgeCallBackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, node v, edgeweight w, edgeid eid):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(u, v, w, eid)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef cppclass TraversalNodeCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(u)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef cppclass TraversalNodeDistCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, count dist):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(u, dist)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef cppclass TraversalNodePairCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, node v):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(u, v)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef extern from "<networkit/graph/BFS.hpp>" namespace "NetworKit::Traversal":

	void BFSfrom[InputIt, Callback](_Graph G, InputIt first, InputIt last, Callback c) nogil except +
	void BFSEdgesFrom[Callback](_Graph G, node source, Callback c) nogil except +

cdef extern from "<networkit/graph/DFS.hpp>" namespace "NetworKit::Traversal":
	void DFSfrom[Callback](_Graph G, node source, Callback c) nogil except +
	void DFSEdgesFrom[Callback](_Graph G, node source, Callback c) nogil except +

cdef class Traversal:

	@staticmethod
	def BFSfrom(Graph graph, start, object callback):
		"""
		Iterate over nodes in breadth-first search order starting from the given node(s).

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.
		start : node/list
			Single node or list of nodes from where the BFS will start.
		callback : Function
			Takes either one (node) or two (node, distance) input parameters.
		"""

		cdef TraversalNodeDistCallbackWrapper *wrapper
		cdef vector[node] sources

		try:
			wrapper = new TraversalNodeDistCallbackWrapper(callback)
			try:
				sources = <vector[node]?>start
			except TypeError:
				sources = [<node?>start]
			BFSfrom[vector[node].iterator, TraversalNodeDistCallbackWrapper](graph._this, sources.begin(),sources.end(), dereference(wrapper))
		finally:
			del wrapper

	@staticmethod
	def BFSEdgesFrom(Graph graph, node start, object callback):
		"""
		Iterate over edges in breadth-first search order starting from the given node(s).

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.
		start : node/list
			Single node or list of nodes from where the BFS will start.
		callback : Function
			Takes four input parameters: (u, v, edgeweight, edgeid)
		"""

		cdef TraversalEdgeCallBackWrapper *wrapper

		try:
			wrapper = new TraversalEdgeCallBackWrapper(callback)
			BFSEdgesFrom[TraversalEdgeCallBackWrapper](graph._this, start, dereference(wrapper))
		finally:
			del wrapper

	@staticmethod
	def DFSfrom(Graph graph, node start, object callback):
		"""
		Iterate over nodes in depth-first search order starting from the given node(s).

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.
		start : node
			Source node from where the DFS will start.
		callback : Function
			Takes a node as input parameter.
		"""
		cdef TraversalNodeCallbackWrapper *wrapper
		try:
			wrapper = new TraversalNodeCallbackWrapper(callback)
			DFSfrom[TraversalNodeCallbackWrapper](graph._this, start, dereference(wrapper))
		finally:
			del wrapper

	@staticmethod
	def DFSEdgesFrom(Graph graph, node start, object callback):
		"""
		Iterate over edges in depth-first search order starting from the given node(s).

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.
		start : node
			Source node from where the DFS will start.
		callback : Function
			Takes four input parameters: (u, v, edgeweight, edgeid)
		"""
		cdef TraversalEdgeCallBackWrapper *wrapper
		try:
			wrapper = new TraversalEdgeCallBackWrapper(callback)
			DFSEdgesFrom[TraversalEdgeCallBackWrapper](graph._this, start, dereference(wrapper))
		finally:
			del wrapper
