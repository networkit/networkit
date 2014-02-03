# cython: language_level=3

#includes


# type imports
from libc.stdint cimport uint64_t
from libc.stdint cimport int64_t


# the C++ standard library
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.map cimport map
from libcpp.string cimport string
from unordered_set cimport unordered_set
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
	#void _configureLogging "Aux::configureLogging" (string loglevel)
	string _getLogLevel "Aux::Log::getLogLevel" () except +
	void _setLogLevel "Aux::Log::setLogLevel" (string loglevel) except +
	void _setPrintLocation "Aux::Log::Settings::setPrintLocation" (bool) except +
	
#def configureLogging(loglevel="ERROR"):
	#""" Set the loglevel of the LOG4CXX module"""
	#_configureLogging(stdstring(loglevel))

def currentLogLevel():
	""" Get the current log level"""
	return pystring(_getLogLevel());

def setLogLevel(loglevel):
	""" Set the current loglevel"""
	_setLogLevel(stdstring(loglevel))

def setPrintLocation(flag):
	""" Switch locations in log statements on or off"""
	_setPrintLocation(flag)

cdef extern from "../src/auxiliary/Parallelism.h" namespace "Aux":
	void _setNumberOfThreads "Aux::setNumberOfThreads" (int)
	int _getCurrentNumberOfThreads "Aux::getCurrentNumberOfThreads" ()
	int _getMaxNumberOfThreads "Aux::getMaxNumberOfThreads" ()
	void _enableNestedParallelism "Aux::enableNestedParallelism" ()

def setNumberOfThreads(nThreads):
	""" Set the number of OpenMP threads """
	_setNumberOfThreads(nThreads)

def getCurrentNumberOfThreads():
	""" Get the number of currently running threads"""
	return _getCurrentNumberOfThreads()

def getMaxNumberOfThreads():
	""" Get the maximum number of available threads"""
	return _getMaxNumberOfThreads()

def enableNestedParallelism():
	""" Enable nested parallelism for OpenMP"""
	_enableNestedParallelism()

# Class definitions


## Module: engineering

# TODO: timer

## Module: graph

cdef extern from "../src/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(count) except +
		count numberOfNodes() except +
		count numberOfEdges() except +
		count degree(node u) except +
		node addNode() except +
		void removeNode(node u) except +
		void addEdge(node u, node v, edgeweight w) except +
		void removeEdge(node u, node v) except +
		bool hasEdge(node u, node v) except +
		edgeweight weight(node u, node v) except +
		vector[node] nodes() except +
		vector[pair[node, node]] edges() except +
		vector[node] neighbors(node u) except +
		void markAsWeighted() except +
		bool isMarkedAsWeighted() except +
		string toString() except +
		string getName() except +
		edgeweight totalEdgeWeight() except +
		

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

	def degree(self, u):
		return self._this.degree(u)
	
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
		
	def neighbors(self, u):
		return self._this.neighbors(u)
	
	def markAsWeighted(self):
		self._this.markAsWeighted()
	
	def isMarkedAsWeighted(self):
		return self._this.isMarkedAsWeighted()

	def toString(self):
		return self._this.toString()

	def getName(self):
		return self._this.getName()

	def totalEdgeWeight(self):
		return self._this.totalEdgeWeight()


cdef class Graph2:
	"""An undirected, optionally weighted graph"""
	cdef _Graph* _this
	
	def __cinit__(self, n=0):
		self._this = new _Graph(n)
		
	# any _thisect which appears as a return type needs to implement setThis
	cdef setThis(self, _Graph* other):
		del self._this
		self._this = other
		return self

	def __dealloc__(self):
		del self._this
	
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



cdef extern from "../src/graph/Subgraph.h":
	cdef cppclass _Subgraph "NetworKit::Subgraph":
		_Subgraph() except +
		_Graph fromNodes(_Graph G, unordered_set[node] nodes)

cdef class Subgraph:
	""" Methods for creating subgraphs"""
	cdef _Subgraph _this

	def fromNodes(self, Graph G, nodes): #unordered_set[node]
		cdef unordered_set[node] nnodes
		for node in nodes:
			nnodes.insert(node);
		return Graph().setThis(self._this.fromNodes(G._this, nnodes))

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


cdef extern from "../src/generators/PubWebGenerator.h":
	cdef cppclass _PubWebGenerator "NetworKit::PubWebGenerator":
		_PubWebGenerator(count numNodes, count numberOfDenseAreas, float neighborhoodRadius, count maxNumberOfNeighbors) except +
		_Graph generate() except +

cdef class PubWebGenerator:
	"""
	 Generates a static graph that resembles an assumed geometric distribution of nodes in
	 a P2P network. The basic structure is to distribute points randomly in the unit torus
	 and to connect vertices close to each other (at most @a neighRad distance and none of
	 them already has @a maxNeigh neighbors). The distribution is chosen to get some areas with
	 high density and others with low density. There are @a numDenseAreas dense areas, which can
	 overlap. Each area is circular, has a certain position and radius and number of points.
	 These values are strored in @a denseAreaXYR and @a numPerArea, respectively.

	 Used and described in more detail in J. Gehweiler, H. Meyerhenke: A Distributed
	 Diffusive Heuristic for Clustering a Virtual P2P Supercomputer. In Proc. 7th High-Performance
	 Grid Computing Workshop (HPGC'10), in conjunction with 24th IEEE Internatl. Parallel and
	 Distributed Processing Symposium (IPDPS'10), IEEE, 2010.

	 Reasonable parameters for constructor:
	 - numNodes: up to a few thousand (possibly more if visualization is not desired and quadratic
	   time complexity has been resolved)
	 - numberOfDenseAreas: [10, 50]
	 - neighborhoodRadius: [0.1, 0.35]
	 - maxNumberOfNeighbors: [4, 40]
	"""
	cdef _PubWebGenerator* _this

	def __cinit__(self, numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors):
		self._this = new _PubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	def generate(self):
		return Graph(0).setThis(self._this.generate())


cdef extern from "../src/generators/ErdosRenyiGenerator.h":
	cdef cppclass _ErdosRenyiGenerator "NetworKit::ErdosRenyiGenerator":
		_ErdosRenyiGenerator(count nNodes, double prob) except +
		_Graph generate() except +

cdef class ErdosRenyiGenerator:
	"""
	  Creates random graphs in the G(n,p) model.
	  The generation follows Vladimir Batagelj and Ulrik Brandes: "Efficient
	  generation of large random networks", Phys Rev E 71, 036113 (2005).
	 
	 Parameters:
	  - nNodes Number of nodes n in the graph.
	  - prob Probability of existence for each edge p.
	"""

	cdef _ErdosRenyiGenerator* _this

	def __cinit__(self, nNodes, prob):
		self._this = new _ErdosRenyiGenerator(nNodes, prob)

	def generate(self):
		return Graph(0).setThis(self._this.generate())


cdef extern from "../src/generators/ChungLuGenerator.h":
	cdef cppclass _ChungLuGenerator "NetworKit::ChungLuGenerator":
		# TODO: revert to count when cython issue fixed
		_ChungLuGenerator(vector[unsigned int] degreeSequence) except +
		_Graph generate() except +

cdef class ChungLuGenerator:
	"""
		Given an arbitrary degree sequence, the Chung-Lu generative model
		will produce a random graph with the same expected degree sequence. 

 		see Aiello, Chung, Lu: A Random Graph Model for Massive Graphs
	"""

	cdef _ChungLuGenerator* _this

	def __cinit__(self, degreeSequence):
		cdef vector
		self._this = new _ChungLuGenerator(degreeSequence)

	def generate(self):
		return Graph(0).setThis(self._this.generate())


cdef extern from "../src/generators/HavelHakimiGenerator.h":
	cdef cppclass _HavelHakimiGenerator "NetworKit::HavelHakimiGenerator":
		# TODO: revert to count when cython issue fixed
		_HavelHakimiGenerator(vector[unsigned int] degreeSequence, bool skipTest) except +
		_Graph generate() except +
		bool isRealizable() except +
		bool getRealizable() except +

cdef class HavelHakimiGenerator:
	"""
	 Havel-Hakimi algorithm for generating a graph according to a given degree sequence.
	 The sequence, if it is realizable, is reconstructed exactly. The resulting graph usually has a high clustering coefficient. Construction runs in linear time O(m). However, the test if a sequence is realizable is quadratic in the sequence length.
	"""

	cdef _HavelHakimiGenerator* _this


	def __cinit__(self, degreeSequence, skipTest=False):
		"""
			@param[in] sequence Degree sequence to realize. Must be non-increasing.
	   		@param[in] skipTest If true, the test if the sequence is realizable is skipped.
	              Default value is false. Set ONLY to true if you are certain that the
	              sequence is realizable
		"""
		self._this = new _HavelHakimiGenerator(degreeSequence, skipTest)

	def isRealizable(self):
		return self._this.isRealizable()

	def getRealizable(self):
		return self._this.getRealizable();

	def generate(self):
		return Graph(0).setThis(self._this.generate())




# Module: graphio

cdef extern from "../src/io/METISGraphReader.h":
	cdef cppclass _METISGraphReader "NetworKit::METISGraphReader":
		_METISGraphReader() except +
		_Graph read(string path) except +
		_Graph* readToHeap(string path) except +

cdef class METISGraphReader:
	""" Reads the METIS adjacency file format [1]
		[1]: http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
	"""
	cdef _METISGraphReader _this
	
	def read(self, path):
		pathbytes = path.encode("utf-8") # string needs to be converted to bytes, which are coerced to std::string
		return Graph(0).setThis(self._this.read(pathbytes))

	def readToHeap(self, path):
		return Graph2(0).setThis(self._this.readToHeap(path.encode("utf-8")))


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

cdef extern from "../src/io/METISGraphWriter.h":
	cdef cppclass _METISGraphWriter "NetworKit::METISGraphWriter":
		_METISGraphWriter() except +
		void write(_Graph G, string path) except +


cdef class METISGraphWriter:
	""" Writes graphs in the METIS format"""
	cdef _METISGraphWriter _this
	
	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		self._this.write(G._this, stdstring(path)) 
	

cdef extern from "../src/io/DotGraphWriter.h":
	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter":
		_DotGraphWriter() except +
		void write(_Graph G, string path) except +


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
		_Graph read(string path) except +
		void write(_Graph G, string path) except +

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
		count numberOfClusters() except +
		cluster clusterOf(node) except +
		float getImbalance() except +
		vector[count] clusterSizes() except +
		map[index, count] clusterSizeMap() except +
		vector[node] getMembers(cluster C) except +

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

	def clusterSizeMap(self):
		return self._this.clusterSizeMap()

	def getMembers(self, C):
		return self._this.getMembers(C)

	def clusterOf(self, v):
		return self._this.clusterOf(v)


cdef extern from "../src/clustering/Coverage.h":
	cdef cppclass _Coverage "NetworKit::Coverage":
		_Coverage() except +
		double getQuality(_Clustering _zeta, _Graph _G) except +

cdef class Coverage:
	""" Coverage is the fraction of intra-community edges """
	cdef _Coverage _this
	
	def getQuality(self, Clustering zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)


cdef extern from "../src/clustering/Modularity.h":
	cdef cppclass _Modularity "NetworKit::Modularity":
		_Modularity() except +
		double getQuality(_Clustering _zeta, _Graph _G) except +
		

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
		_PLP(count updateThreshold) except +
		_Clustering run(_Graph _G) except +
		count numberOfIterations() except +
		string toString() except +


cdef class PLP(Clusterer):
	""" Parallel label propagation for community detection: 
		Moderate solution quality, very short time to solution.
	 """
	cdef _PLP _this

	def __cinit__(self, updateThreshold=None):
		if updateThreshold is None:
			self._this = _PLP()
		else:
			self._this = _PLP(updateThreshold)

	
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
		
		
cdef extern from "../src/community/PLM2.h":
	cdef cppclass _PLM2 "NetworKit::PLM2":
		_PLM2() except +
		_PLM2(bool refine, double gamma, string par, count maxIter) except +
		string toString() except +
		_Clustering run(_Graph G) except +


cdef class PLM2(Clusterer):
	""" MultiLevel Parallel LocalMover - the Louvain method, optionally extended to
		a full multi-level algorithm with refinement"""
		
	cdef _PLM2 _this
	
	def __cinit__(self, refine=True, gamma=1.0, par="balanced", maxIter=32):
		self._this = _PLM2(refine, gamma, stdstring(par), maxIter)
		
	def toString(self):
		return self._this.toString().decode("utf-8")
		
	def run(self, Graph G not None):
		return Clustering().setThis(self._this.run(G._this))


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

cdef extern from "../src/community/CommunityGraph.h":
	cdef cppclass _CommunityGraph "NetworKit::CommunityGraph":
		void run(_Graph G, _Clustering zeta) except +
		_Graph getGraph() except +
		map[index, node] getCommunityToNodeMap() except +
		map[node, index] getNodeToCommunityMap() except +

cdef class CommunityGraph:
	""" Coarsen a graph according to communities """
	cdef _CommunityGraph _this

	def run(self, Graph G, Clustering zeta):
		self._this.run(G._this, zeta._this)

	def getGraph(self):
		return Graph().setThis(self._this.getGraph())

	def getCommunityToNodeMap(self):
		return self._this.getCommunityToNodeMap()

	def getNodeToCommunityMap(self):
		return self._this.getNodeToCommunityMap()

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
		void run(_Graph G) except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		vector[node] getComponent(index component) except +
		map[index, count] getComponentSizes() except +


cdef class ConnectedComponents:
	""" Determines the connected components and associated values for
		an undirected graph.
	"""
	cdef _ConnectedComponents _this

	def run(self, Graph G):
		self._this.run(G._this)

	def numberOfComponents(self):
		return self._this.numberOfComponents()

	def componentOfNode(self, v):
		return self._this.componentOfNode(v)

	def getComponent(self, componentIndex):
		return self._this.getComponent(componentIndex)
	
	def getComponentSizes(self):
		return self._this.getComponentSizes()



# Module: dynamic

cdef extern from "../src/dynamics/GraphEvent.h":
	enum _Type "NetworKit::GraphEvent::Type": 		
		NODE_ADDITION,
		NODE_REMOVAL,
		EDGE_ADDITION,
		EDGE_REMOVAL,
		EDGE_WEIGHT_UPDATE,
		TIME_STEP

cdef extern from "../src/dynamics/GraphEvent.h":
	cdef cppclass _GraphEvent "NetworKit::GraphEvent":
		node u, v
		edgeweight w
		_Type type
		_GraphEvent() except +
		_GraphEvent(_Type type, node u, node v, edgeweight w) except +
		string toString() except +
		
cdef class GraphEvent:
	cdef _GraphEvent _this

	property type:
		def __get__(self): 
			return self._this.type
		def __set__(self, t):
			self._this.type = t

	property u:
		def __get__(self): 
			return self._this.u
		def __set__(self, u):
			self._this.u = u

	property v:
		def __get__(self): 
			return self._this.v
		def __set__(self, v):
			self._this.type = v

	property w:
		def __get__(self): 
			return self._this.w
		def __set__(self, w):
			self._this.type = w			
	
	def __cinit__(self, type, u, v, w):
		self._this = _GraphEvent(type, u, v, w)

	def toString(self):
		return self._this.toString().decode("utf-8")

	def __repr__(self):
		return self.toString()


cdef extern from "../src/dynamics/DGSStreamParser.h":
	cdef cppclass _DGSStreamParser "NetworKit::DGSStreamParser":
		_DGSStreamParser(string path, bool mapped, node baseIndex) except +
		vector[_GraphEvent] getStream() except +

cdef class DGSStreamParser:
	cdef _DGSStreamParser* _this

	def __cinit__(self, path, mapped=True, baseIndex=0):
		self._this = new _DGSStreamParser(stdstring(path), mapped, baseIndex)

	def getStream(self):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.getStream()]


cdef extern from "../src/dynamics/DGSWriter.h":
	cdef cppclass _DGSWriter "NetworKit::DGSWriter":
		void write(vector[_GraphEvent] stream, string path) except +


cdef class DGSWriter:
	cdef _DGSWriter* _this

	def __cinit__(self):
		self._this = new _DGSWriter()

	def write(self, stream, path):
		cdef vector[_GraphEvent] _stream
		for ev in stream:
			_stream.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.write(_stream, stdstring(path))


cdef extern from "../src/dcd2/DynamicCommunityDetection.h":
	cdef cppclass _DynamicCommunityDetection "NetworKit::DynamicCommunityDetection":
		_DynamicCommunityDetection(string inputPath, string algoName, string updateStrategy, count interval, vector[string] recordSettings) except +
		void run() except +
		vector[double] getTimeline(string key) except +
		vector[pair[count, count]] getGraphSizeTimeline() except +
		vector[pair[_Graph, _Clustering]] getResultTimeline() except +

cdef class DynamicCommunityDetection:
	cdef _DynamicCommunityDetection* _this

	def __cinit__(self, inputPath, algoName, updateStrategy, interval, recordSettings):
		self._this = new _DynamicCommunityDetection(stdstring(inputPath), stdstring(algoName), stdstring(updateStrategy), interval, [stdstring(key) for key in recordSettings])

	def run(self):
		self._this.run()

	def getTimeline(self, key):
		return self._this.getTimeline(stdstring(key))

	def getGraphSizeTimeline(self):
		return self._this.getGraphSizeTimeline()

	def getResultTimeline(self):
		timeline = []
		for pair in self._this.getResultTimeline():
			_G = pair.first
			_zeta = pair.second
			timeline.append((Graph().setThis(_G), Clustering().setThis(_zeta)))
		return timeline
			


cdef extern from "../src/generators/DynamicPathGenerator.h":
	cdef cppclass _DynamicPathGenerator "NetworKit::DynamicPathGenerator":
		_DynamicPathGenerator() except +
		vector[_GraphEvent] generate(count nSteps) except +


cdef class DynamicPathGenerator:
	cdef _DynamicPathGenerator* _this

	def __cinit__(self):
		self._this = new _DynamicPathGenerator()

	def generate(self, nSteps):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]





cdef extern from "../src/generators/DynamicPubWebGenerator.h":
	cdef cppclass _DynamicPubWebGenerator "NetworKit::DynamicPubWebGenerator":
		_DynamicPubWebGenerator(count numNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors) except +
		vector[_GraphEvent] generate(count nSteps) except +
		_Graph getGraph() except +


cdef class DynamicPubWebGenerator:
	cdef _DynamicPubWebGenerator* _this

	def __cinit__(self, numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors):
		self._this = new _DynamicPubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	def generate(self, nSteps):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]

	def getGraph(self):
		return Graph().setThis(self._this.getGraph())


cdef extern from "../src/generators/ForestFireGenerator.h":
	cdef cppclass _ForestFireGenerator "NetworKit::ForestFireGenerator":
		_ForestFireGenerator(double p) except +
		vector[_GraphEvent] generate(count nSteps) except +
		_Graph getGraph() except +


cdef class ForestFireGenerator:
	cdef _ForestFireGenerator* _this

	def __cinit__(self, p):
		self._this = new _ForestFireGenerator(p)

	def generate(self, nSteps):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]



cdef extern from "../src/dynamics/GraphUpdater.h":
	cdef cppclass _GraphUpdater "NetworKit::GraphUpdater":
		_GraphUpdater(_Graph G) except +
		void update(vector[_GraphEvent] stream) except +
		vector[pair[count, count]] getSizeTimeline() except +

cdef class GraphUpdater:
	cdef _GraphUpdater* _this

	def __cinit__(self, Graph G):
		self._this = new _GraphUpdater(G._this)

	def update(self, stream):
		cdef vector[_GraphEvent] _stream
		for ev in stream:
			_stream.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.update(_stream)


