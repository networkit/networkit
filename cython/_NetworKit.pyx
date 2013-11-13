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




# Cython helper functions

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes

def pystring(stdstring):
	""" convert a std::string (= python byte string) to a normal Python string"""
	return stdstring.decode("utf-8")


# Function definitions

cdef extern from "../src/auxiliary/Log.h" namespace "Aux":
	void _configureLogging "Aux::configureLogging" (string loglevel)
	string _currentLogLevel "Aux::currentLogLevel" ()
	void _setLoglevel "Aux::setLoglevel" (string loglevel)
	
def configureLogging(loglevel="ERROR"):
	""" Set the loglevel of the LOG4CXX module"""
	_configureLogging(stdstring(loglevel))

def currentLogLevel():
	""" Get the current log level"""
	return pystring(_currentLogLevel());

def setLoglevel(loglevel):
	""" Set the current loglevel"""
	_setLoglevel(stdstring(loglevel))


cdef extern from "../src/auxiliary/Parallelism.h" namespace "Aux":
	void _setNumberOfThreads "Aux::setNumberOfThreads" (int)
	int _getCurrentNumberOfThreads "Aux::getCurrentNumberOfThreads" ()
	int _getMaxNumberOfThreads "Aux::getMaxNumberOfThreads" ()

def setNumberOfThreads(nThreads):
	""" Set the number of OpenMP threads """
	_setNumberOfThreads(nThreads)

def getCurrentNumberOfThreads():
	""" Get the number of currently running threads"""
	return _getCurrentNumberOfThreads()

def getMaxNumberOfThreads():
	""" Get the maximum number of available threads"""
	return _getMaxNumberOfThreads()

# Class definitions


## Module: engineering

# TODO: timer

## Module: graph

cdef extern from "../src/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(count) except +
		count numberOfNodes()
		count numberOfEdges()
		node addNode()
		void removeNode(node u)
		void addEdge(node u, node v, edgeweight w)
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
	"""An undirected, optionally weighted graph"""
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


cdef extern from "../src/graph/BFS.h":
	cdef cppclass _BFS "NetworKit::BFS":
		_BFS() except +
		vector[count] run(_Graph G, node source)

cdef class BFS:
	""" Simple breadth-first search"""
	cdef _BFS _this

	def run(self, Graph G not None, source):
		"""
		Breadth-first search from source.
		return Vector of unweighted distances from node @a source, i.e. the
	 		length (number of edges) of the shortest path from source to any other vertex.
		"""
		return self._this.run(G._this, source)
	

cdef extern from "../src/graph/Dijkstra.h":
	cdef cppclass _Dijkstra "NetworKit::Dijkstra":
		_Dijkstra() except +
		vector[edgeweight] run(_Graph G, node source)

cdef class Dijkstra:
	""" Dijkstra's SSSP algorithm.
	 	Returns list of weighted distances from node source, i.e. the
	    length of the shortest path from @a source to any other node."""
	cdef _Dijkstra _this

	def run(self, Graph G not None, source):
		return self._this.run(G._this, source)


cdef extern from "../src/independentset/Luby.h":
	cdef cppclass _Luby "NetworKit::Luby":
		_Luby() except +
		vector[bool] run(_Graph G)
		string toString()

# FIXME: check correctness
cdef class Luby:
	""" Luby's parallel maximal independent set algorithm"""
	cdef _Luby _this

	def run(self, Graph G not None):
		return self._this.run(G._this)

	def toString(self):
		return self._this.toString().decode("utf-8")


# Module: generators
	
cdef extern from "../src/graph/GraphGenerator.h":
	cdef cppclass _GraphGenerator "NetworKit::GraphGenerator":
		_GraphGenerator() except +
		_Graph makeRandomGraph(count n, double p)


cdef class GraphGenerator:
	""" Provides several functions for graph generation"""
	cdef _GraphGenerator _this
	
	def __cinit__(self):
		self._this = _GraphGenerator()
		
	
	def makeRandomGraph(self, n, p):
		cdef _Graph _G = self._this.makeRandomGraph(n, p)
		return Graph(0).setThis(_G)

cdef extern from "../src/generators/BarabasiAlbertGenerator.h":
	cdef cppclass _BarabasiAlbertGenerator "NetworKit::BarabasiAlbertGenerator":
		_BarabasiAlbertGenerator() except +
		_BarabasiAlbertGenerator(count k, count nMax, count n0) except +
		_Graph generate()

cdef class BarabasiAlbertGenerator:
	""" Generates a scale-free graph according to the Barabasi-Albert model """
	cdef _BarabasiAlbertGenerator _this

	def __cinit__(self, k, nMax, n0):
		self._this = _BarabasiAlbertGenerator(k, nMax, n0)

	def generate(self):
		return Graph().setThis(self._this.generate());



# Module: graphio

cdef extern from "../src/io/METISGraphReader.h":
	cdef cppclass _METISGraphReader "NetworKit::METISGraphReader":
		_METISGraphReader() except +
		_Graph read(string path)

cdef class METISGraphReader:
	""" Reads the METIS adjacency file format [1]
		[1]: http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
	"""
	cdef _METISGraphReader _this
	
	def read(self, path):
		pathbytes = path.encode("utf-8") # string needs to be converted to bytes, which are coerced to std::string
		cdef _Graph _G = self._this.read(pathbytes)
		return Graph(0).setThis(_G)


cdef extern from "../src/io/FastMETISGraphReader.h":
	cdef cppclass _FastMETISGraphReader "NetworKit::FastMETISGraphReader":
		_FastMETISGraphReader() except +
		_Graph read(string path) except +

cdef class FastMETISGraphReader:
	""" A faster but currently experimental implementation of a reader for
		the METIS format.
	"""
	cdef _FastMETISGraphReader _this
	
	def read(self, path):
		pathbytes = path.encode("utf-8") # string needs to be converted to bytes, which are coerced to std::string
		cdef _Graph _G = self._this.read(pathbytes)
		return Graph(0).setThis(_G)
	
	

cdef extern from "../src/io/DotGraphWriter.h":
	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter":
		_DotGraphWriter() except +
		void write(_Graph G, string path)


cdef class DotGraphWriter:
	""" Writes graphs in the .dot/GraphViz format"""
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
	""" Reads and writes graphs in various edge list formats"""

	cdef _EdgeListIO _this

	def __cinit__(self, separator, firstNode):
		cdef char sep = stdstring(separator)[0]
		self._this = _EdgeListIO(sep, firstNode)

	def read(self, path):
		return Graph().setThis(self._this.read(stdstring(path)))

	def write(self, Graph G not None, path):
		self._this.write(G._this, stdstring(path))



cdef extern from "../src/io/LineFileReader.h":
	cdef cppclass _LineFileReader "NetworKit::LineFileReader":
		_LineFileReader() except +
		vector[string] read(string path)


cdef class LineFileReader:
	""" Reads a file and puts each line in a list of strings """
	cdef _LineFileReader _this

	def read(self, path):
		return self._this.read(stdstring(path))


cdef extern from "../src/io/SNAPGraphWriter.h":
	cdef cppclass _SNAPGraphWriter "NetworKit::SNAPGraphWriter":
		_SNAPGraphWriter() except +
		void write(_Graph G, string path) 

cdef class SNAPGraphWriter:
	""" Writes graphs in a format suitable for the Georgia Tech SNAP software [1]
		[1]: http://snap-graph.sourceforge.net/
	"""
	cdef _SNAPGraphWriter _this

	def write(self, Graph G, path):
		self._this.write(G._this, stdstring(path))


cdef extern from "../src/io/ClusteringReader.h":
	cdef cppclass _ClusteringReader "NetworKit::ClusteringReader":
		_ClusteringReader() except +
		_Clustering read(string path)


cdef class ClusteringReader:
	""" Reads a partition from a file.
		File format: line i contains subset id of element i.
	 """
	cdef _ClusteringReader _this

	def read(self, path):
		return Clustering().setThis(self._this.read(stdstring(path)))


cdef extern from "../src/io/ClusteringWriter.h":
	cdef cppclass _ClusteringWriter "NetworKit::ClusteringWriter":
		_ClusteringWriter() except +
		void write(_Clustering, string path)


cdef class ClusteringWriter:
	""" Writes a partition to a file.
		File format: line i contains subset id of element i.
	 """
	cdef _ClusteringWriter _this

	def write(self, Clustering zeta, path):
		self._this.write(zeta._this, stdstring(path))


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



# Module: community


cdef extern from "../src/clustering/Clustering.h":
	cdef cppclass _Clustering "NetworKit::Clustering":
		_Clustering() except +
		count numberOfClusters()
		cluster clusterOf(node)
		float getImbalance()
		vector[count] clusterSizes()
		vector[node] getMembers(cluster C)

cdef class Clustering:
	"""
		Represents a partition of the node set into disjoint subsets (clusters).
	"""
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

	def clusterOf(self, v):
		return self._this.clusterOf(v)


cdef extern from "../src/clustering/Coverage.h":
	cdef cppclass _Coverage "NetworKit::Coverage":
		_Coverage() except +
		double getQuality(_Clustering _zeta, _Graph _G)

cdef class Coverage:
	""" Coverage is the fraction of intra-community edges """
	cdef _Coverage _this
	
	def getQuality(self, Clustering zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)


cdef extern from "../src/clustering/Modularity.h":
	cdef cppclass _Modularity "NetworKit::Modularity":
		_Modularity() except +
		double getQuality(_Clustering _zeta, _Graph _G)
		

cdef class Modularity:
	"""
		Modularity [1] is a quality index for community detection.
		[1]: http://en.wikipedia.org/wiki/Modularity_%28networks%29
	"""
	cdef _Modularity _this
	
	def getQuality(self, Clustering zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)


cdef class Clusterer:
	""" Abstract base class for static community detection algorithms"""
	pass

cdef extern from "../src/community/PLP.h":
	cdef cppclass _PLP "NetworKit::PLP":
		_PLP() except +
		_Clustering run(_Graph _G)
		count numberOfIterations()
		string toString()


cdef class PLP(Clusterer):
	""" Parallel label propagation for community detection: 
		Moderate solution quality, very short time to solution.
	 """
	cdef _PLP _this
	
	def run(self, Graph G not None):
		return Clustering().setThis(self._this.run(G._this))

	def numberOfIterations(self):
		return self._this.numberOfIterations()

	def toString(self):
		return self._this.toString().decode("utf-8")


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
	

cdef extern from "../src/community/PLM.h":
	cdef cppclass _PLM "NetworKit::PLM":
		_PLM() except +
		_PLM(string par, double gamma)
		_Clustering run(_Graph _G)
		string toString()
		
cdef class PLM(Clusterer):
	""" Parallel Louvain method for community detection: 
	High solution quality, moderate time to solution. """

	cdef _PLM _this
	
	def __cinit__(self, par="balanced", gamma=1.0):
		self._this = _PLM(stdstring(par), gamma)
	
	def run(self, Graph G not None):
		return Clustering().setThis(self._this.run(G._this))

	def toString(self):
		return self._this.toString().decode("utf-8")


# FIXME: PLM2 
# FIXME: CNM

# # PLM2

# cdef extern from "../src/community/PLM2.h":
# 	cdef cppclass _PLM2 "NetworKit::PLM2":
# 		_PLM2() except +
# 		_PLM2(string par, double gamma)
# 		_Clustering run(_Graph _G)
# 		string toString()

# cdef class PLM2:
# 	cdef _PLM2 _this

# 	def __cinit__(self, par="balanced", gamma=1.0):
# 		self._this = _PLM2(stdstring(par), gamma)

# 	def run(self, Graph G):
# 		return Clustering().setThis(self._this.run(G._this))

# 	def toString(self):
# 		return self._this.toString().decode("utf-8")


# cdef extern from "../src/community/CNM.h":
# 	cdef cppclass _CNM "NetworKit::CNM":
# 		_CNM() except +
# 		_Clustering run(_Graph _G)
# 		string toString()

# cdef class CNM:
# 	cdef _CNM _this

# 	def run(self, Graph G):
# 		return Clustering().setThis(self._this.run(G._this))

# 	def toString(self):
# 		return self._this.toString().decode("utf-8")


cdef class DissimilarityMeasure:
	""" Abstract base class for partition/community dissimilarity measures"""
	pass


cdef extern from "../src/clustering/NodeStructuralRandMeasure.h":
	cdef cppclass _NodeStructuralRandMeasure "NetworKit::NodeStructuralRandMeasure":
		_NodeStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Clustering first, _Clustering second)

cdef class NodeStructuralRandMeasure(DissimilarityMeasure):
	""" The node-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering all pairs of nodes.
	"""
	cdef _NodeStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Clustering first, Clustering second):
		return self._this.getDissimilarity(G._this, first._this, second._this)


cdef extern from "../src/clustering/GraphStructuralRandMeasure.h":
	cdef cppclass _GraphStructuralRandMeasure "NetworKit::GraphStructuralRandMeasure":
		_GraphStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Clustering first, _Clustering second)

cdef class GraphStructuralRandMeasure(DissimilarityMeasure):
	""" The graph-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering connected pairs of nodes.
	"""
	cdef _GraphStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Clustering first, Clustering second):
		return self._this.getDissimilarity(G._this, first._this, second._this)


cdef extern from "../src/community/EPP.h":
	cdef cppclass _EPP "NetworKit::EPP":
		_Clustering run(_Graph G)
		string toString()

cdef class EPP(Clusterer):
	""" EPP - Ensemble Preprocessing """
	cdef _EPP _this

	def run(self, Graph G):
		return Clustering().setThis(self._this.run(G._this))

	def toString(self):
		return self._this.toString()

	cdef setThis(self, _EPP other):
		self._this = other
		return self


cdef extern from "../src/community/EPPFactory.h":
	cdef cppclass _EPPFactory "NetworKit::EPPFactory":
		_EPP make(count ensembleSize, string baseAlgorithm, string finalAlgorithm)

cdef class EPPFactory:
	""" This class makes instaces of the EPP community detection algorithm"""
	cdef _EPPFactory _this
	
	def make(self, ensembleSize, baseAlgorithm="PLP", finalAlgorithm="PLM"):
		return EPP().setThis(self._this.make(ensembleSize, stdstring(baseAlgorithm), stdstring(finalAlgorithm)))



# Module: properties

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
	""" Collects various functions for basic graph properties """

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


cdef extern from "../src/properties/ConnectedComponents.h":
	cdef cppclass _ConnectedComponents "NetworKit::ConnectedComponents":
		_ConnectedComponents() except +
		void run(_Graph G)
		count numberOfComponents()
		count sizeOfComponent(index component)
		count componentOfNode(node query)
		vector[node] getComponent(index component)
		vector[count] getComponentSizes()


cdef class ConnectedComponents:
	""" Determines the connected components and associated values for
		an undirected graph.
	"""
	cdef _ConnectedComponents _this

	def run(self, Graph G):
		self._this.run(G._this)

	def numberOfComponents(self):
		return self._this.numberOfComponents()

	def sizeOfComponent(self, componentIndex):
		return self._this.sizeOfComponent(componentIndex)

	def componentOfNode(self, v):
		return self._this.componentOfNode(v)

	def getComponent(self, componentIndex):
		return self._this.getComponent(componentIndex)
	
	def getComponentSizes(self):
		return self._this.getComponentSizes()






