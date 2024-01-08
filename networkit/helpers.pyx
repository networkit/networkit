# distutils: language=c++

from cython.operator cimport preincrement, dereference

from libcpp.vector cimport vector
from libcpp.utility cimport pair

import numpy as np

cdef extern from "<networkit/auxiliary/Parallel.hpp>" namespace "Aux::Parallel":

	void sort[Iter](Iter begin, Iter end) nogil
	void sort[Iter, Comp](Iter begin, Iter end, Comp compare) nogil

class ReducedFunctionalityWarning(Warning):
	"Reduced networkit functionality warning"
	pass

def ranked(sample):
	"""
	ranked(sample)

	Given a list of numbers, this function computes the rank of each value
	and returns a list of ranks where result[i] is the rank of
	the i-th element in the given sample.
	Currently used in networkit.profiling.stat.

	Parameters
	----------
	sample : list(float)
		Input list of values.

	Returns
	-------
	list(float)
		Ranking of input values.
	"""
	cdef vector[pair[double, count]] helper = vector[pair[double, count]](len(sample))
	cdef vector[double] result = vector[double](len(sample), 0)
	for i in range(len(sample)):
		helper[i] = <pair[double, count]?>(sample[i], i)
	sort(helper.begin(), helper.end())
	cdef double value = helper[0].first
	cdef double summ = 0.
	cdef count length = 0
	for i in range(len(sample)):
		if value == helper[i].first:
			summ += (i+1)
			length += 1
		else:
			summ /= length
			for j in range(length):
				result[helper[i-j-1].second] = summ
			value = helper[i].first
			summ = i+1
			length = 1
	summ /= length
	for j in range(length):
		result[helper[len(sample)-j-1].second] = summ
	return result

# Cython helper functions
def stdstring(pystring):
	""" 
	stdstring(pystring)

	Convert a Python string to a bytes object which is automatically coerced to std::string.

	Parameters
	----------
	pystring : str
		Input python string.

	Returns
	-------
	stdstring
		Python bytes string.
	"""
	return pystring.encode("utf-8")

def pystring(stdstring):
	""" 
	pystring(stdstring)

	Convert a std::string (= python byte string) to a normal Python string.

	Parameters
	----------
	stdstring : str
		Input python byte string.

	Returns
	-------
	pystring
		Python string.
	"""
	return stdstring.decode("utf-8")


cdef get_dtype(element_t* _=NULL):
	"""
	Infer the numpy dtype from a fused cython type.
	"""
	if element_t is float:
		return np.float32
	elif element_t is double:
		return np.float64
	elif element_t is int:
		return np.int32
	elif element_t is long:
		return np.int64
	else:
		raise TypeError


cdef asarray_1d(vector[element_t]* vec):
	"""
	Convert a vector to a one-dimensional array.

	Parameters
	----------
	vec : vector[element_t]
		Vector representing a one-dimensional array.

	Returns
	-------
	np.ndarray
	"""
	cdef:
		element_t[:] values
		element_t* target

	values = np.empty(vec.size(), get_dtype[element_t]())
	target = &values[0]
	iterator = vec.begin()
	while iterator != vec.end():
		target[0] = dereference(iterator)
		preincrement(iterator)
		preincrement(target)

	return values


cdef asarray_2d(vector[vector[element_t]]* nested):
	"""
	Convert nested vectors to a two-dimensional array.

	Parameters
	----------
	nested : vector[vector[element_t]]
		Nested vectors representing a two-dimensional array.

	Returns
	-------
	np.ndarray
	"""
	cdef:
		element_t[:] values
		element_t* target
		int num_rows, num_cols

	# Return an empty matrix if there are no rows.
	num_rows = nested.size()
	if num_rows == 0:
		return np.empty((num_rows, num_rows), get_dtype[element_t]())

	# Allocate memory.
	row_iterator = nested.begin()
	num_cols = dereference(row_iterator).size()
	values = np.empty(num_rows * num_cols, get_dtype[element_t]())

	# Populate the memory.
	target = &values[0]
	while row_iterator != nested.end():
		row = dereference(row_iterator)
		iterator = row.begin()
		while iterator != row.end():
			target[0] = dereference(iterator)
			preincrement(iterator)
			preincrement(target)
		preincrement(row_iterator)
	return np.reshape(values, (num_rows, num_cols))
