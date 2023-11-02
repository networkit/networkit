#distutils: language=c++

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph

ctypedef uint64_t count

cdef extern from "<networkit/subgraphs/GraphletsCounter.hpp>":
	cdef cppclass _GraphletsCounter "NetworKit::GraphletsCounter"(_Algorithm):
		_GraphletsCounter(_Graph G, unsigned K)  # noexcept

		const vector[count]& getGraphletsCounts() const  # noexcept

cdef class GraphletsCounter(Algorithm):
	""" Count the number of occurrences of every graphets
	of a given size in an undirected graph.

	GraphletsCounter(G, K)

	Create a `K`-graphlets counter for Graph `G`.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	K : int
		The graphlet size.
	"""
	cdef Graph _G
	cdef unsigned _K

	def __cinit__(self, Graph G, unsigned K):
		self._G = G
		self._K = K
		self._this = new _GraphletsCounter(G._this, K)

	def getCounts(self):
		"""
		"""
		return (<_GraphletsCounter*>(self._this)).getGraphletsCounts()
