from libcpp cimport bool as bool_t
from libcpp.string cimport string

cdef extern from "<networkit/base/Algorithm.hpp>" namespace "NetworKit":
	cdef cppclass _Algorithm "NetworKit::Algorithm":
		_Algorithm()
		void run() nogil except +
		bool_t hasFinished() except +
		string toString() except +
		bool_t isParallel() except +

cdef class Algorithm:
	cdef _Algorithm *_this
