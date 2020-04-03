from libc.stdint cimport uint64_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool as bool_t

ctypedef uint64_t index
ctypedef index node

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition

cdef extern from "<networkit/centrality/Centrality.hpp>":

	cdef cppclass _Centrality "NetworKit::Centrality"(_Algorithm):
		_Centrality(_Graph, bool_t, bool_t) except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +
		double maximum() except +
		double centralization() except +


cdef class Centrality(Algorithm):
	""" Abstract base class for centrality measures"""

	cdef Graph _G
