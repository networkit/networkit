# type imports
from libc.stdint cimport int64_t

ctypedef int64_t count

# Cython class definitions

cdef extern from "../src/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(int) except +
		count numberOfNodes()
		

cdef extern from "../src/graph/GraphGenerator.h":
	cdef cppclass _GraphGenerator "NetworKit::GraphGenerator":
		_GraphGenerator() except +
		_Graph makeRandomGraph(count n, double p)

		
cdef class Graph:
	cdef _Graph obj
	
	
	def __cinit__(self, n):
		self.obj = _Graph(n)
		
# 	def __dealloc__(self):
# 		del self.obj
	
	# any object which appears as a return type needs to implement setThis
	cdef setThis(self, _Graph other):
		#del self.obj
		self.obj = other
	
	def numberOfNodes(self):
		return self.obj.numberOfNodes()


cdef class GraphGenerator:
	cdef _GraphGenerator obj
	
	def __cinit__(self):
		self.obj = _GraphGenerator()
		
	
# 	def __dealloc__(self):
# 		del self.obj
	
	def makeRandomGraph(self, n, p):
		cdef _Graph _G = self.obj.makeRandomGraph(n, p)
		G = Graph(0)
		G.setThis(_G)
		return G

#cdef _Graph* test = &_Graph(0)
