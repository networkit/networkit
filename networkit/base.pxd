from libcpp cimport bool as bool_t

cdef extern from "<networkit/base/Algorithm.hpp>" namespace "NetworKit":
	cdef cppclass _Algorithm "NetworKit::Algorithm":
		_Algorithm()
		void run() except + nogil
		bool_t hasFinished() except +

cdef class _CythonParentClass:
	cdef _Algorithm *_this

cdef class Algorithm(_CythonParentClass):
	pass
