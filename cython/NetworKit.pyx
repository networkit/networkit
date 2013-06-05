# type imports
from libc.stdint cimport int64_t

from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string

ctypedef int64_t count
ctypedef int64_t index
ctypedef index node


# python module imports
import networkx as nx


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
		vector[node] nodes()
		vector[pair[node, node]] edges()
		

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
		
	def nodes(self):
		return self.obj.nodes()
	
	def edges(self):
		return self.obj.edges()
	
	
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
		return self


cdef extern from "../src/community/LabelPropagation.h":
	cdef cppclass _LabelPropagation "NetworKit::LabelPropagation":
		_LabelPropagation() except +
		_Clustering run(_Graph _G)

cdef class LabelPropagation:
	cdef _LabelPropagation obj
	
	def run(self, Graph G not None):
		return Clustering().setObj(self.obj.run(G.obj))


cdef extern from "../src/io/DotGraphWriter.h":
	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter":
		_DotGraphWriter() except +
		void write(_Graph G, string path)


cdef class DotGraphWriter:
	cdef _DotGraphWriter obj
	
	def write(self, Graph G not None, path):
		pathbytes = path.encode("utf-8") # string needs to be converted to bytes, which are coerced to std::string
		self.obj.write(G.obj, pathbytes)
	

cdef extern from "../src/viz/ForceDirected.h":
	cdef cppclass _ForceDirected "NetworKit::ForceDirected":
		_ForceDirected() except +
		void draw(_Graph _G)
	
cdef class ForceDirected:
	cdef _ForceDirected obj
	
	def draw(self, Graph G not None):
		pass
	


# TODO: initialize log4cxx

#--------- NetworKit functions ----------------#

def readGraph(path):
	"""	Automatically detect input format and read graph"""
	# TODO: detect file format by looking at the file content
	if path.endswith(".graph"):
		reader = METISGraphReader()
	else:
		raise Exception("unknown graph file format")
	G = reader.read(path)
	return G




def nx2nk(nxG):
	""" Convert a networkx.Graph to a NetworKit.Graph """
	
	n = nxG.number_of_nodes()
	cdef Graph nkG = Graph(n)
	
	for (u, v) in nxG.edges():
		nkG.addEdge(u, v)
	
	return nkG

def nk2nx(nkG):
	""" Convert a NetworKit.Graph to a networkx.Graph """
	nxG = nx.Graph()
	for (u, v) in nkG.edges():
		nxG.add_edge(u, v)
	
	return nxG
	
	