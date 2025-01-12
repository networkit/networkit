from cython.operator cimport dereference

from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

from .graph cimport _Graph
from .structures cimport _Partition, _Cover
from .structures cimport count, index, node, edgeweight

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Graph move( _Graph t ) nogil # specialized declaration as general declaration disables template argument deduction and doesn't work
	_Partition move( _Partition t) nogil
	_Cover move(_Cover t) nogil
	vector[double] move(vector[double])
	vector[bool_t] move(vector[bool_t])
	vector[count] move(vector[count])
	pair[_Graph, vector[node]] move(pair[_Graph, vector[node]]) nogil
	vector[pair[pair[node, node], double]] move(vector[pair[pair[node, node], double]]) nogil
	vector[pair[node, node]] move(vector[pair[node, node]]) nogil

ctypedef fused element_t:
	edgeweight
	node
	double

cdef asarray_1d(vector[element_t]* vec)
cdef asarray_2d(vector[vector[element_t]]* nested)

cdef inline maybe_asarray_1d(vector[element_t]* vec, asarray):
	return asarray_1d[element_t](vec) if asarray else dereference(vec)

cdef inline maybe_asarray_2d(vector[vector[element_t]]* nested, asarray):
	return asarray_2d[element_t](nested) if asarray else dereference(nested)
