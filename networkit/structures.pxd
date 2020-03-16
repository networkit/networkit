from libc.stdint cimport uint64_t
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string

ctypedef uint64_t count
ctypedef uint64_t index

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Partition move( _Partition t) nogil
	_Cover move(_Cover t) nogil

cdef extern from "<networkit/structures/Cover.hpp>":

	cdef cppclass _Cover "NetworKit::Cover":
		_Cover() except +
		_Cover(_Partition p) except +
		_Cover(count n) except +
		set[index] subsetsOf(index e) except +
		index extend() except +
		void remove(index e) except +
		void addToSubset(index s, index e) except +
		void removeFromSubset(index s, index e) except +
		void moveToSubset(index s, index e) except +
		void toSingleton(index e) except +
		void allToSingletons() except +
		void mergeSubsets(index s, index t) except +
		void setUpperBound(index upper) except +
		index upperBound() except +
		index lowerBound() except +
		bool_t contains(index e) except +
		bool_t inSameSubset(index e1, index e2) except +
		vector[count] subsetSizes() except +
		map[index, count] subsetSizeMap() except +
		set[index] getMembers(const index s) except +
		count numberOfElements() except +
		count numberOfSubsets() except +
		set[index] getSubsetIds() except +

cdef class Cover:
	cdef _Cover _this
	cdef setThis(self, _Cover& other)

cdef extern from "<networkit/structures/Partition.hpp>":

	cdef cppclass _Partition "NetworKit::Partition":
		_Partition() except +
		_Partition(index) except +
		_Partition(_Partition) except +
		_Partition(vector[index]) except +
		index subsetOf(index e) except +
		index extend() except +
		void remove(index e) except +
		void addToSubset(index s, index e) except +
		void moveToSubset(index s, index e) except +
		void toSingleton(index e) except +
		void allToSingletons() except +
		index mergeSubsets(index s, index t) except +
		void setUpperBound(index upper) except +
		index upperBound() except +
		index lowerBound() except +
		void compact(bool_t useTurbo) except +
		bool_t contains(index e) except +
		bool_t inSameSubset(index e1, index e2) except +
		vector[count] subsetSizes() except +
		map[index, count] subsetSizeMap() except +
		set[index] getMembers(const index s) except +
		count numberOfElements() except +
		count numberOfSubsets() except +
		vector[index] getVector() except +
		void setName(string name) except +
		string getName() except +
		set[index] getSubsetIds() except +
		index operator[](index) except +

cdef class Partition:
	cdef _Partition _this
	cdef setThis(self, _Partition& other)
