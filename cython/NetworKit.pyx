# type imports
from libc.stdint cimport int64_t

from libcpp.vector cimport vector
from libcpp.string cimport string

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
		
	# any object which appears as a return type needs to implement setObj
	cdef setObj(self, _Graph other):
		#del self.obj
		self.obj = other
		return self
	
	def numberOfNodes(self):
		return self.obj.numberOfNodes()


cdef class GraphGenerator:
	cdef _GraphGenerator obj
	
	def __cinit__(self):
		self.obj = _GraphGenerator()
		
	
	def makeRandomGraph(self, n, p):
		cdef _Graph _G = self.obj.makeRandomGraph(n, p)
		return Graph(0).setObj(_G)



cdef extern from "../src/io/METISGraphReader.h":
	cdef cppclass _METISGraphReader "NetworKit::METISGraphReader":
		_METISGraphReader() except +
		_Graph read(string path)

cdef class METISGraphReader:
	cdef _METISGraphReader obj
	
	def read(self, path):
		pathbytes = path.encode("utf-8") # string needs to be converted to bytes, which are coerced to std::string
		cdef _Graph _G = self.obj.read(pathbytes)
		return Graph(0).setObj(_G)
