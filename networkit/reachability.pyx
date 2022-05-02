# distutils: language=c++

from cython.operator import dereference, preincrement

from libcpp.vector cimport vector
from libcpp cimport bool as bool_t
from libcpp.string cimport string

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .helpers import stdstring
from .structures cimport count, index, node

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

none = _none

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef extern from "<networkit/reachability/ReachableNodes.hpp>":

	cdef cppclass _ReachableNodes "NetworKit::ReachableNodes"(_Algorithm):
		_ReachableNodes(_Graph G, bool_t exact) except +
		count numberOfReachableNodes(node u) except +
		count numberOfReachableNodesLB(node u) except +
		count numberOfReachableNodesUB(node u) except +
		bool_t exact

cdef class ReachableNodes(Algorithm):
	"""
	ReachableNodes(G, exact)	

	Determines or estimates the number of reachable nodes from each node in the graph.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	exact : bool
		Whether or not to compute the number of reachable nodes exactly. Only
		used for directed graphs, on undirected graphs the number of reachable
		nodes from every node can be computed in linear time.
	"""

	cdef Graph _G

	def __cinit__(self, Graph G, exact = True):
		self._G = G
		self._this = new _ReachableNodes(G._this, exact)

	def numberOfReachableNodes(self, node u):
		"""
		numberOfReachableNodes(u)

		Returns the number of reachable nodes from the given node 'u'. Only available if 'exact' is true.

		Parameters
		----------
		u : int
			A node.

		Returns
		-------
		int
			The number of nodes reachable from 'u'.
		"""
		return (<_ReachableNodes*>(self._this)).numberOfReachableNodes(u)

	def numberOfReachableNodesLB(self, node u):
		"""
		numberOfReachableNodesLB(u)

		Returns a lower bound of the number of reachable nodes from the given node 'u'.

		Parameters
		----------
		u : int
			A node.

		Returns
		-------
		int
			Lower bound of number of nodes reachable from 'u'.
		"""
		return (<_ReachableNodes*>(self._this)).numberOfReachableNodesLB(u)

	def numberOfReachableNodesUB(self, node u):
		"""
		numberOfReachableNodesUB(u)

		Returns an upper bound of the number of reachable nodes from the given node 'u'.

		Parameters
		----------
		u : int
			A node.

		Returns
		-------
		int
			Upper bound of number of nodes reachable from 'u'.
		"""
		return (<_ReachableNodes*>(self._this)).numberOfReachableNodesUB(u)

	property exact:
		def __get__(self):
			return (<_ReachableNodes*>(self._this)).exact
		def __set__(self, bool_t exact):
			""" Use a different edge direction. """
			(<_ReachableNodes*>(self._this)).exact = exact

cdef cppclass PathCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(vector[node] path):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(path)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef extern from "<networkit/reachability/AllSimplePaths.hpp>":

	cdef cppclass _AllSimplePaths "NetworKit::AllSimplePaths"(_Algorithm):
		_AllSimplePaths(_Graph G, node source, node target, count cutoff) except +
		void run() nogil except +
		count numberOfSimplePaths() except +
		vector[vector[node]] getAllSimplePaths() except +
		void forAllSimplePaths[Callback](Callback c) except +

cdef class AllSimplePaths(Algorithm):
	"""
	AllSimplePaths(G, source, target, cutoff=None)

	Algorithm to compute all existing simple paths from a source node to a
	target node. The maximum length of the paths can be fixed through 'cutoff'.

	Note
	----
	CAUTION: This algorithm could take a lot of time on large networks (many
	edges), especially if the cutoff value is high or not specified.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	source : int
		The source node.
	target : int
		The target node.
	cutoff : int, optional
		The maximum length of the simple paths. Default: None

	"""

	cdef Graph _G

	def __cinit__(self, Graph G, source, target, cutoff=none):
		self._G = G
		self._this = new _AllSimplePaths(G._this, source, target, cutoff)

	def numberOfSimplePaths(self):
		"""
		numberOfSimplePaths()

		Returns the number of simple paths.

		Returns
		-------
		int
			The number of simple paths.
		"""
		return (<_AllSimplePaths*>(self._this)).numberOfSimplePaths()

	def getAllSimplePaths(self):
		"""
		getAllSimplePaths()

		Returns all the simple paths from source to target.

		Returns
		-------
		list(list(int))
			A list containing all simple paths (each represented by a list of nodes).
		"""
		return (<_AllSimplePaths*>(self._this)).getAllSimplePaths()

	def forAllSimplePaths(self, object callback):
		""" 
		forAllSimplePaths(callback)		

		More efficient path iterator. Iterates over all the simple paths.

		Parameters
		-----------
		callback : object
			Any callable object that takes the parameter path
		"""
		cdef PathCallbackWrapper* wrapper
		try:
			wrapper = new PathCallbackWrapper(callback)
			(<_AllSimplePaths*>(self._this)).forAllSimplePaths[PathCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper
