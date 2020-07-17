from libcpp.string cimport string
from libc.stdint cimport uint64_t
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node

from .graph cimport _Graph
from .structures cimport _Partition, _Cover
from .matching cimport _Matching

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Graph move( _Graph t ) nogil # specialized declaration as general declaration disables template argument deduction and doesn't work
	_Partition move( _Partition t) nogil
	_Cover move(_Cover t) nogil
	_Matching move(_Matching) nogil
	vector[double] move(vector[double])
	vector[bool_t] move(vector[bool_t])
	vector[count] move(vector[count])
	pair[_Graph, vector[node]] move(pair[_Graph, vector[node]]) nogil
	vector[pair[pair[node, node], double]] move(vector[pair[pair[node, node], double]]) nogil
	vector[pair[node, node]] move(vector[pair[node, node]]) nogil

