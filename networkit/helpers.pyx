# distutils: language=c++

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

cdef extern from "<networkit/auxiliary/Parallel.hpp>" namespace "Aux::Parallel":

	void sort[Iter](Iter begin, Iter end) nogil
	void sort[Iter, Comp](Iter begin, Iter end, Comp compare) nogil

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

def sorted(sample):
	"""	
	sorted(sample)	

	Returns a sorted list of given numbers.

	Note
	----
	DEPRECATED. Use :code:`sorted()` function provided by Python.

	Parameters
	----------
	sample : list(float)
		(Unsorted) input list of values.

	Returns
	-------
	list(float)
		Sorted list of values.
	"""
	cdef vector[double] result = <vector[double]?>sample
	sort(result.begin(),result.end())
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

