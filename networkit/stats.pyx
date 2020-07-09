# distutils: language=c++

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

ctypedef uint64_t count

cdef extern from "<networkit/auxiliary/Parallel.hpp>" namespace "Aux::Parallel":

	void sort[Iter](Iter begin, Iter end) nogil
	void sort[Iter, Comp](Iter begin, Iter end, Comp compare) nogil
""" Missing statistics functions"""

# stats

def gini(values):
	"""
	Computes the Gini coefficient for the distribution given as a list of values.
	"""
	sorted_list = sorted(values)
	height, area = 0, 0
	for value in sorted_list:
		height += value
		area += height - value / 2.
	fair_area = height * len(values) / 2
	return (fair_area - area) / fair_area



