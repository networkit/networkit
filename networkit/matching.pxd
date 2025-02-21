from libcpp.vector cimport vector
from libcpp.unordered_set cimport unordered_set
from libcpp cimport bool as bool_t
from libcpp.string cimport string

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition, count, index, node, edgeweight

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Matching move(_Matching) nogil
	_BMatching move(_BMatching) nogil

cdef extern from "<networkit/matching/Matching.hpp>":

	cdef cppclass _Matching "NetworKit::Matching":
		_Matching() except +
		_Matching(count) except +
		void match(node, node) except +
		void unmatch(node, node) except +
		bool_t isMatched(node) except +
		bool_t areMatched(node, node) except +
		bool_t isProper(_Graph) except +
		count size(_Graph) except +
		index mate(node) except +
		edgeweight weight(_Graph) except +
		_Partition toPartition(_Graph) except +
		vector[node] getVector() except +

cdef class Matching:
	cdef _Matching _this
	cdef setThis(self,  _Matching& other)

cdef extern from "<networkit/matching/BMatching.hpp>":

	cdef cppclass _BMatching "NetworKit::BMatching":
		_BMatching() except +
		_BMatching(_Graph, vector[count]) except +
		bool_t isProper() except +
		void match(node, node) except +
		void unmatch(node, node) except +
		bool_t isUnmatched(node) except +
		bool_t areMatched(node, node) except +
		count size() except +
		edgeweight weight() except +
		vector[unordered_set[node]] &getMatches() except +
		vector[count] getB() except +
		void reset() except +

cdef class BMatching:
	cdef _BMatching _this
	cdef Graph _G
	cdef setThis(self,  _BMatching& other)