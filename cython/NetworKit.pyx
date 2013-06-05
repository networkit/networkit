# type imports
from libc.stdint cimport int64_t

from libcpp.vector cimport vector
from libcpp.string cimport string

ctypedef int64_t count
ctypedef int64_t index
ctypedef index node

# Cython class definitions

cdef extern from "../src/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(int) except +
		count numberOfNodes()
		count numberOfEdges()
		node addNode()
		void removeNode(node u)
		void addEdge(node u, node v)
		void removeEdge(node u, node v)
		

cdef class Graph:
	cdef _Graph obj
	
	def __cinit__(self, n):
		self.obj = _Graph(n)
		
	# any object which appears as a return type needs to implement setObj
	cdef setObj(self, _Graph other):
		#del self.obj
		self.obj = other
		return self
	
	cdef _Graph getObj(self):
		return self.obj
	
	def numberOfNodes(self):
		return self.obj.numberOfNodes()
	
	def numberOfEdges(self):
		return self.obj.numberOfEdges()
	
	def addNode(self):
		return self.obj.addNode()
	
	def removeNode(self, u):
		self.obj.removeNode(u)
		
	def addEdge(self, u, v):
		self.obj.addEdge(u, v)
		
	def removeEdge(self, u, v):
		self.obj.removeEdge(u, v)
	
	
cdef extern from "../src/graph/GraphGenerator.h":
	cdef cppclass _GraphGenerator "NetworKit::GraphGenerator":
		_GraphGenerator() except +
		_Graph makeRandomGraph(count n, double p)


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
	

cdef extern from "../src/clustering/Clustering.h":
	cdef cppclass _Clustering "NetworKit::Clustering":
		_Clustering() except +

cdef class Clustering:
	cdef _Clustering obj
	
	cdef setObj(self, _Clustering other):
		self.obj = other


cdef extern from "../src/community/LabelPropagation.h":
	cdef cppclass _LabelPropagation "NetworKit::LabelPropagation":
		_LabelPropagation() except +
		_Clustering run(_Graph _G)

cdef class LabelPropagation:
	cdef _LabelPropagation obj
	
	def run(self, Graph G not None):
		self.obj.run(G.obj)


# TODO: initialize log4cxx


def readGraph(path):
	# TODO: detect file format by looking at the file content
	if path.endswith(".graph"):
		reader = METISGraphReader()
	else:
		raise Exception("unknown graph file format")
	G = reader.read(path)
	return G