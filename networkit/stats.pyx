# distutils: language=c++

from libcpp.vector cimport vector
from libcpp.utility cimport pair

from .structures cimport count

cdef extern from "<networkit/auxiliary/Parallel.hpp>" namespace "Aux::Parallel":

	void sort[Iter](Iter begin, Iter end) nogil
	void sort[Iter, Comp](Iter begin, Iter end, Comp compare) nogil
""" Missing statistics functions"""

# stats

def gini(values):
	"""
	gini(values)

	Computes the Gini coefficient for the distribution given as a list of values.

	Parameters
	----------
	values : list(float)
		Input list containing values.

	Returns
	-------
	float
		Gini coefficient. Empty or all-zero input yields 0.0.
	"""
	if len(values) == 0:
		return 0.0
	sorted_list = sorted(values)
	height, area = 0, 0
	for value in sorted_list:
		height += value
		area += height - value / 2.
	fair_area = height * len(values) / 2
	if fair_area == 0:
		return 0.0
	return (fair_area - area) / fair_area
