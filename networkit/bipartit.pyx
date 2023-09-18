# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.unordered_set cimport unordered_set

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition, count, index, node

cdef extern from "<networkit/bipartit/Bipartit.hpp>":

	cdef cppclass _Bipartit "NetworKit::Bipartit"(_Algorithm):
		_Bipartit(_Graph G) except +
		bool_t isBipartit() except +
		_Partition getPartition() except +
		vector[node] getOddCycle() except +

cdef class Bipartit(Algorithm):
	"""
	Bipartit()

	Determines if a graph is bipartit.
	If the graph is bipartit a bipartition is provided.
	Otherwise an odd cycle is provided.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""

	def __cinit__(self, Graph G):
		self._this = new _Bipartit(G._this)

	def isBipartit(self):
		"""
		isBipartit()

		Returns if the graph is bipartit.

		Returns
		-------
		bool
			true iff the graph is bipartit
		"""
		return (<_Bipartit*>(self._this)).isBipartit()

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
		return Partition().setThis((<_Bipartit*>(self._this)).getPartition())

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
		return (<_Bipartit*>(self._this)).getOddCycle()