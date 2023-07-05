from libcpp.vector cimport vector

from .dynamics cimport _GraphEvent, GraphEvent

cdef extern from "<networkit/base/DynAlgorithm.hpp>":
	cdef cppclass _DynAlgorithm "NetworKit::DynAlgorithm":
		void update(_GraphEvent) nogil except +
		void updateBatch(vector[_GraphEvent]) nogil except +
