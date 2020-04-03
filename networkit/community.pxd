from libc.stdint cimport uint64_t

from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.map cimport map

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition 

cdef extern from "<algorithm>" namespace "std":
	pair[_Graph, vector[node]] move(pair[_Graph, vector[node]]) nogil

cdef extern from "<networkit/community/CommunityDetectionAlgorithm.hpp>":

	cdef cppclass _CommunityDetectionAlgorithm "NetworKit::CommunityDetectionAlgorithm"(_Algorithm):
		_CommunityDetectionAlgorithm(const _Graph &_G)
		_Partition getPartition() except +

cdef class CommunityDetector(Algorithm):
	cdef Graph _G
