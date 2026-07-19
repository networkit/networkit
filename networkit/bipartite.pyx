# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.unordered_set cimport unordered_set

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition, count, index, node

cdef extern from "<networkit/bipartite/Bipartite.hpp>":

	cdef cppclass _Bipartite "NetworKit::Bipartite"(_Algorithm):
		_Bipartite(_Graph G) except +
		bool_t isBipartite() except +
		_Partition getPartition() except +
		vector[node] getOddCycle() except +

cdef class Bipartite(Algorithm):
	"""
	Bipartite()

	Determines if a graph is bipartite.
	If the graph is bipartite a bipartition is provided.
	Otherwise an odd cycle is provided.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""

	def __cinit__(self, Graph G):
		self._this = new _Bipartite(G._this)

	def isBipartite(self):
		"""
		isBipartite()

		Returns if the graph is bipartite.

		Returns
		-------
		bool
			true iff the graph is bipartite
		"""
		return (<_Bipartite*>(self._this)).isBipartite()

	def getPartition(self):
		"""
		getPartition()

		Returns a bipartition of the graph if one exists.
		Throws error otherwise.

		Returns
		-------
		Partition
			bipartition of the graph (throws error if none exists)
		"""
		return Partition().setThis((<_Bipartite*>(self._this)).getPartition())

	def getOddCycle(self):
		"""
		getOddCycle()

		Returns an odd cycle if one exists.
		Throws error otherwise.

		Returns
		-------
		vector[node]
			A vector of nodes that form an odd cycle when read in order with the forst and last element being connected.
		"""
		return (<_Bipartite*>(self._this)).getOddCycle()