from libcpp.vector cimport vector

from .dynamics cimport _GraphEvent, GraphEvent

cdef extern from "<networkit/base/DynAlgorithm.hpp>":
	cdef cppclass _DynAlgorithm "NetworKit::DynAlgorithm":
		void update(_GraphEvent) except + nogil
		void updateBatch(vector[_GraphEvent]) except + nogil
