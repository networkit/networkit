# type imports
from libc.stdint cimport int64_t

ctypedef int64_t count

# Cython class definitions

cdef extern from "../src/auxiliary/RandomInteger.h":
	cdef cppclass _RandomInteger "Aux::RandomInteger":
		_RandomInteger() except +
		int64_t generate(int64_t lower, int64_t upper)
		

cdef extern from "../src/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph(int) except +
		count numberOfNodes()
		

cdef extern from "../src/graph/GraphGenerator.h":
	cdef cppclass _GraphGenerator "NetworKit::GraphGenerator":
		_GraphGenerator() except +
		_Graph makeRandomGraph(count n, double p)

		
		
# Python wrappers
	
cdef class RandomInteger:

	cdef _RandomInteger *thisptr
	
	def __cinit__(self):
		self.thisptr = new _RandomInteger()
		
	def __dealloc__(self):
		del self.thisptr
	
	def generate(self, lower, upper):
		return self.thisptr.generate(lower, upper)

cdef class Graph:

	cdef _Graph *thisptr
	
	def __cinit__(self, n):
		self.thisptr = new _Graph(n)
		
	def __dealloc__(self):
		del self.thisptr
		
	cdef setInstance(self, _Graph *ptr):
		del self.thisptr
		self.thisptr = ptr
	
	def numberOfNodes(self):
		return self.thisptr.numberOfNodes()


cdef class GraphGenerator:

	cdef _GraphGenerator *thisptr
	
	def __cinit__(self):
		self.thisptr = new _GraphGenerator()
		
	
	def __dealloc__(self):
		del self.thisptr
	
	def makeRandomGraph(self, n, p):
		G = Graph(0)
		#cdef _Graph* _G = &(self.thisptr.makeRandomGraph(n, p))
		#G.setInstance(_G)
		return G

#cdef _Graph* test = &_Graph(0)