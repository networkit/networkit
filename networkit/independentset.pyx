# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.string cimport string

from .graph cimport _Graph, Graph

# TODO: expose all methods
cdef extern from "<networkit/independentset/Luby.hpp>":

	cdef cppclass _Luby "NetworKit::Luby":
		_Luby() except +
		vector[bool_t] run(_Graph G) except +
		string toString()


# FIXME: check correctness
cdef class Luby:
	""" Luby's parallel maximal independent set algorithm"""
	cdef _Luby _this

	def run(self, Graph G not None):
		""" Returns a bool vector of length n where vec[v] is True iff v is in the independent sets.
		
		Parameters:
		-----------
		G : networkit.Graph
			The graph.
		
		Returns:
		--------
		vector
			A bool vector of length n.
		"""
		return self._this.run(G._this)
		# TODO: return self

	def toString(self):
		""" Get string representation of the algorithm.
		
		Returns:
		--------
		string
			The string representation of the algorithm.
		"""
		return self._this.toString().decode("utf-8")

