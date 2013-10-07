# cython: language_level=3

#includes


# type imports
from libc.stdint cimport uint64_t
from libc.stdint cimport int64_t


# the C++ standard library
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
# FIXME: from libcpp.unordered_map import unordered_map

# NetworKit typedefs
ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node
ctypedef index cluster
ctypedef double edgeweight




# helper functions

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes


# Cython class definitions

cdef extern from "../src/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(count) except +
		count numberOfNodes()
		count numberOfEdges()
		node addNode()
		void removeNode(node u)
		void addEdge(node u, node v, edgeweight w)
		# TODO: optional weight argument
		void removeEdge(node u, node v)
		bool hasEdge(node u, node v)
		edgeweight weight(node u, node v)
		vector[node] nodes()
		vector[pair[node, node]] edges()
		void markAsWeighted()
		bool isMarkedAsWeighted()
		string toString()
		string getName()
		

cdef class Graph:
	cdef _Graph _this
	
	def __cinit__(self, n=None):
		if n is not None:
			self._this = _Graph(n)
		
	# any _thisect which appears as a return type needs to implement setThis
	cdef setThis(self, _Graph other):
		#del self._this
		self._this = other
		return self
	
	def numberOfNodes(self):
		return self._this.numberOfNodes()
	
	def numberOfEdges(self):
		return self._this.numberOfEdges()
	
	def addNode(self):
		return self._this.addNode()
	
	def removeNode(self, u):
		self._this.removeNode(u)
		
	def addEdge(self, u, v, w=1.0):
		self._this.addEdge(u, v, w)
		
	def removeEdge(self, u, v):
		self._this.removeEdge(u, v)
		
	def hasEdge(self, u, v):
		return self._this.hasEdge(u, v)
		
	def weight(self, u, v):
		return self._this.weight(u, v)
		
	def nodes(self):
		return self._this.nodes()
	
	def edges(self):
		return self._this.edges()
	
	def markAsWeighted(self):
		self._this.markAsWeighted()
	
	def isMarkedAsWeighted(self):
		return self._this.isMarkedAsWeighted()

	def toString(self):
		return self._this.toString()

	def getName(self):
		return self._this.getName()
	
	
cdef extern from "../src/graph/GraphGenerator.h":
	cdef cppclass _GraphGenerator "NetworKit::GraphGenerator":
		_GraphGenerator() except +
		_Graph makeRandomGraph(count n, double p)


cdef class GraphGenerator:
	cdef _GraphGenerator _this
	
	def __cinit__(self):
		self._this = _GraphGenerator()
		
	
	def makeRandomGraph(self, n, p):
		cdef _Graph _G = self._this.makeRandomGraph(n, p)
		return Graph(0).setThis(_G)



cdef extern from "../src/io/METISGraphReader.h":
	cdef cppclass _METISGraphReader "NetworKit::METISGraphReader":
		_METISGraphReader() except +
		_Graph read(string path)

cdef class METISGraphReader:
	cdef _METISGraphReader _this
	
	def read(self, path):
		pathbytes = path.encode("utf-8") # string needs to be converted to bytes, which are coerced to std::string
		cdef _Graph _G = self._this.read(pathbytes)
		return Graph(0).setThis(_G)
	

cdef extern from "../src/io/DotGraphWriter.h":
	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter":
		_DotGraphWriter() except +
		void write(_Graph G, string path)


cdef class DotGraphWriter:
	cdef _DotGraphWriter _this
	
	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		self._this.write(G._this, stdstring(path))
		

cdef extern from "../src/io/EdgeListIO.h":
	cdef cppclass _EdgeListIO "NetworKit::EdgeListIO":
		_EdgeListIO() except +
		_EdgeListIO(char separator, node firstNode) except +
		_Graph read(string path)
		void write(_Graph G, string path)

cdef class EdgeListIO:

	cdef _EdgeListIO _this

	def __cinit__(self, separator, firstNode):
		cdef char sep = stdstring(separator)[0]
		self._this = _EdgeListIO(sep, firstNode)

	def read(self, path):
		return Graph().setThis(self._this.read(stdstring(path)))

	def write(self, Graph G not None, path):
		self._this.write(G._this, stdstring(path))

# cdef extern from "../src/io/EdgeListIO.h":
# 	cdef cppclass _EdgeListIO "NetworKit::EdgeListIO":
# 		_EdgeListIO() except +
# 		_EdgeListIO(char separator, node firstNode) except +
# 		_Graph read(string path)
# 		void write(_Graph G, string path)

# cdef class EdgeListIO:

# 	cdef _EdgeListIO _this
	
# 	def __cinit__(self, string separator, node firstNode):
# 		cdef char sep = separator[0]
# 		self._this = _EdgeListIO(sep, firstNode)
		
# 	def read(self, path):
# 		return Graph().setThis(self._this.read(stdstring(path)))
	
# 	def write(self, Graph G not None, path):
# 		 # string needs to be converted to bytes, which are coerced to std::string
# 		self._this.write(G._this, stdstring(path))




cdef extern from "../src/clustering/Clustering.h":
	cdef cppclass _Clustering "NetworKit::Clustering":
		_Clustering() except +
		count numberOfClusters()
		float getImbalance()
		vector[count] clusterSizes()
		vector[node] getMembers(cluster C)

cdef class Clustering:
	cdef _Clustering _this
	
	cdef setThis(self, _Clustering other):
		self._this = other
		return self

	def numberOfClusters(self):
		return self._this.numberOfClusters()

	def getImbalance(self):
		return self._this.getImbalance()

	def clusterSizes(self):
		return self._this.clusterSizes()

	def getMembers(self, C):
		return self._this.getMembers(C)



cdef class Clusterer:
	""" Abstract base class for static community detection algorithms"""
	def run(self, Graph G not None):
		raise NotImplementedError("abstract method")

cdef extern from "../src/community/LabelPropagation.h":
	cdef cppclass _LabelPropagation "NetworKit::LabelPropagation":
		_LabelPropagation() except +
		_Clustering run(_Graph _G)
		count numberOfIterations()


cdef class LabelPropagation(Clusterer):
	cdef _LabelPropagation _this
	
	def run(self, Graph G not None):
		return Clustering().setThis(self._this.run(G._this))

	def numberOfIterations(self):
		return self._this.numberOfIterations()


cdef extern from "../src/community/LPDegreeOrdered.h":
	cdef cppclass _LPDegreeOrdered "NetworKit::LPDegreeOrdered":
		_LPDegreeOrdered() except +
		_Clustering run(_Graph _G)
		count numberOfIterations()

cdef class LPDegreeOrdered(Clusterer):
	cdef _LPDegreeOrdered _this
	
	def run(self, Graph G not None):
		return Clustering().setThis(self._this.run(G._this))

	def numberOfIterations(self):
		return self._this.numberOfIterations()
	
	

cdef extern from "../src/community/Louvain.h":
	cdef cppclass _Louvain "NetworKit::Louvain":
		_Louvain() except +
		_Louvain(string par, double gamma)
		_Clustering run(_Graph _G)
		
cdef class Louvain(Clusterer):
	cdef _Louvain _this
	
	def __cinit__(self, par="balanced", gamma=1.0):
		self._this = _Louvain(stdstring(par), gamma)
	
	def run(self, Graph G not None):
		return Clustering().setThis(self._this.run(G._this))



# this is an example for using static methods
cdef extern from "../src/properties/GraphProperties.h" namespace "NetworKit::GraphProperties":
	# static methods live in the class namespace, so declare them here
	pair[count, count] minMaxDegree(_Graph _G)
	double averageDegree(_Graph _G)
	vector[count] degreeDistribution(_Graph _G)
	vector[double] localClusteringCoefficients(_Graph _G)
	double averageLocalClusteringCoefficient(_Graph _G)
	vector[double] localClusteringCoefficientPerDegree(_Graph _G)
	
	cdef cppclass _GraphProperties "NetworKit::GraphProperties":
		pass

cdef class GraphProperties:

	@staticmethod
	def minMaxDegree(Graph G not None):
		return minMaxDegree(G._this)

	@staticmethod
	def averageDegree(Graph G not None):
		return averageDegree(G._this)

	@staticmethod
	def degreeDistribution(Graph G not None):
		return degreeDistribution(G._this)

	@staticmethod
	def averageLocalClusteringCoefficient(Graph G not None):
		return averageLocalClusteringCoefficient(G._this)


cdef extern from "../src/dynamics/GraphEventHandler.h":
	cdef cppclass _GraphEventHandler "NetworKit::GraphEventHandler":
		_GraphEventHandler() except +
		void onNodeAddition(node u)
		void onNodeRemoval(node u)
		void onEdgeAddition(node u, node v)
		void onEdgeRemoval(node u, node v)
		void onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew)
		void onTimeStep()


cdef extern from "../src/dynamics/GraphEventProxy.h":
	cdef cppclass _GraphEventProxy "NetworKit::GraphEventProxy":
		_GraphEventProxy()	# nullary constructor not valid
		_GraphEventProxy(_Graph _G)
		void registerObserver(_GraphEventHandler* _observer)
		node addNode()
		void removeNode(node u)
		void addEdge(node u, node v)
		# TODO: optional edge weight
		void removeEdge(node u, node v)
		void setWeight(node u, node v, edgeweight w)
		void timeStep()
		
cdef class GraphEventProxy:
	cdef _GraphEventProxy _this
	
	def __cinit__(self, Graph G not None):
		self._this = _GraphEventProxy(G._this)
	# TODO: delegates



cdef extern from "../src/io/DGSReader.h":
	cdef cppclass _DGSReader "NetworKit::DGSReader":
		_DGSReader() except +
		void read(string path, _GraphEventProxy _proxy)
		
		
cdef class DGSReader:
	cdef _DGSReader _this
	
	def read(self, path, GraphEventProxy proxy not None):
		self._this.read(stdstring(path), proxy._this)
		
cdef extern from "../src/clustering/Modularity.h":
	cdef cppclass _Modularity "NetworKit::Modularity":
		_Modularity() except +
		double getQuality(_Clustering _zeta, _Graph _G)
		

cdef class Modularity:
	cdef _Modularity _this
	
	def getQuality(self, Clustering zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)



cdef extern from "../src/io/SNAPGraphWriter.h":
	cdef cppclass _SNAPGraphWriter "NetworKit::SNAPGraphWriter":
		_SNAPGraphWriter() except +
		void write(_Graph G, string path) 

cdef class SNAPGraphWriter:
	cdef _SNAPGraphWriter _this

	def write(self, Graph G, path):
		self._this.write(G._this, stdstring(path))

# cdef extern from "../src/dcd/DynamicLabelPropagation.h":
# 	cdef cppclass _DynamicLabelPropagation "NetworKit::DynamicLabelPropagation":
# 		_DynamicLabelPropagation() except +
# 		_DynamicLabelPropagation(_Graph G, count theta, string strategyName) except +
# 		_Clustering run()
# 		string toString()
# 
# class DynamicCommunityDetector:
# 	pass	
# 		
# cdef class DynamicLabelPropagation:
# 	cdef _DynamicLabelPropagation _this
# 	
# 	def __cinit__(self, Graph G not None, theta, strategyName):
# 		self._this = _DynamicLabelPropagation(G._this, theta, stdstring(strategyName))
# 		
# 	def run(self):
# 		self._this.run()


# FIXME:
# cdef extern from "../src/generators/DynamicBarabasiAlbertGenerator.h":
# 	cdef cppclass _DynamicBarabasiAlbertGenerator "NetworKit::DynamicBarabasiAlbertGenerator":
# 		_DynamicBarabasiAlbertGenerator() except +
# 		_DynamicBarabasiAlbertGenerator(_GraphEventProxy _Gproxy, count k) except +
# 		void initializeGraph()
# 		void generate()
# 		
# cdef class DynamicBarabasiAlbertGenerator:
# 	cdef _DynamicBarabasiAlbertGenerator _this
# 	
# 	def __cinit__(self, GraphEventProxy Gproxy not None, k):
# 		self._this = _DynamicBarabasiAlbertGenerator(Gproxy._this, k)
# 		
# 	def generate(self):
# 		self._this.generate()
		


# under construction

cdef extern from "../src/generators/PubWebGenerator.h":
	cdef cppclass _PubWebGenerator "NetworKit::PubWebGenerator":
		_PubWebGenerator() except +
		_PubWebGenerator(count numNodes, count numberOfDenseAreas, float neighborhoodRadius, count maxNumberOfNeighbors) except +
		_Graph generate()
		
cdef class PubWebGenerator:
	cdef _PubWebGenerator _this
	
	def __cinit__(self, count numNodes, count numberOfDenseAreas, float neighborhoodRadius, count maxNumberOfNeighbors):
		self._this = _PubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)
		
	def generate(self):
		return Graph().setThis(self._this.generate())	

cdef extern from "../src/viz/ForceDirected.h":
	cdef cppclass _ForceDirected "NetworKit::ForceDirected":
		_ForceDirected() except +
		void draw(_Graph _G)
	
cdef class ForceDirected:
	cdef _ForceDirected _this
	
	def draw(self, Graph G not None):
		pass

# Parameters

cdef extern from "../src/base/Parameters.h":
	cdef cppclass _Parameters "NetworKit::Parameters":
		_Parameters() except +
		void setInt(string key, int64_t value)
		void setDouble(string key, double value)
		void setString(key, value)
		void setBool(string key, bool value)
		int64_t getInt(string key)
		double getDouble(string key)
		string getString(string key)
		bool getBool(string key)

# TODO: wrapper for Parameters


# SelectiveSCAN

# cdef extern from "../src/scd/SelectiveSCAN.h":
# 	cdef cppclass _TSelectiveSCAN "NetworKit::TSelectiveSCAN":
# 		_TSelectiveSCAN() except +
# 		_TSelectiveSCAN(_Graph G, _Parameters param, double epsilon, double mu) except +
# 		unordered_map[node, pair[unordered_set[node], int64_t]] run(unordered_set[node] seeds)


# cdef class TSelectiveSCAN:
# 	pass


# LineFileReader

cdef extern from "../src/io/LineFileReader.h":
	cdef cppclass _LineFileReader "NetworKit::LineFileReader":
		_LineFileReader() except +
		vector[string] read(string path)


cdef class LineFileReader:
	cdef _LineFileReader _this

	def read(self, path):
		return self._this.read(stdstring(path))

# TODO: initialize log4cxx

	