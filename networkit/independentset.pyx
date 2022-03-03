# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.string cimport string

from .graph cimport _Graph, Graph

cdef extern from "<networkit/independentset/IndependentSetFinder.hpp>":

	cdef cppclass _IndependentSetFinder "NetworKit::IndependentSetFinder":
		_IndependentSetFinder() except +
		vector[bool_t] run(const _Graph& G) except +
		bool_t isIndependentSet(const vector[bool_t]& nodeSet, const _Graph& G)

cdef class IndependentSetFinder:
	""" 
	IndependentSetFinder()
	
	Abstract base class for independent set algorithms.
	"""
	cdef _IndependentSetFinder* _this

	def __cinit__(self, *args):
		# The construction is handled by the subclasses
		return

	def __dealloc__(self):
		if self._this is not NULL:
			del self._this
			self._this = NULL

	def run(self, Graph G):
		"""
		run(G)

		Returns a list, where an entry indicates whether node i is part of an independent set.

		Parameters
		----------
		G : networkit.Graph
			The graph.

		Returns
		-------
		list(bool)
			A list of bools.
		"""
		return self._this.run(G._this)

	def isIndependentSet(self, nodeSet, Graph G):
		""" 
		isIndependentSet(set, G)
		
		Checks whether a given set is an independent set in graph G.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		nodeSet : list(bool)
			Input list, where each entry represents a node. If entry is True, the node is part of the set; False if not.

		Returns
		-------
		bool

		"""
		return self._this.isIndependentSet(nodeSet, G._this)

cdef extern from "<networkit/independentset/Luby.hpp>":

	cdef cppclass _Luby "NetworKit::Luby"(_IndependentSetFinder):
		_Luby() except +

cdef class Luby(IndependentSetFinder):
	""" 
	Luby()

	Luby's parallel maximal independent set algorithm.
	"""

	def __cinit__(self, *args):
		self._this = new _Luby()
