# cython: language_level=3

#includes

# C++ operators
from cython.operator import dereference

# type imports
from libc.stdint cimport uint64_t
from libc.stdint cimport int64_t


# the C++ standard library
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string
from unordered_set cimport unordered_set
from unordered_map cimport unordered_map

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

cdef extern from "../cpp/auxiliary/Log.h" namespace "Aux":
	#void _configureLogging "Aux::configureLogging" (string loglevel)
	string _getLogLevel "Aux::Log::getLogLevel" () except +
	void _setLogLevel "Aux::Log::setLogLevel" (string loglevel) except +
	void _setPrintLocation "Aux::Log::Settings::setPrintLocation" (bool) except +

def getLogLevel():
	""" Get the current log level"""
	return pystring(_getLogLevel());

def setLogLevel(loglevel):
	""" Set the current loglevel"""
	_setLogLevel(stdstring(loglevel))

def setPrintLocation(flag):
	""" Switch locations in log statements on or off"""
	_setPrintLocation(flag)

cdef extern from "../cpp/auxiliary/Parallelism.h" namespace "Aux":
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

cdef extern from "../cpp/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(count, bool, bool) except +
		void stealFrom(_Graph)
		count numberOfNodes() except +
		count numberOfEdges() except +
		count degree(node u) except +
		count degreeIn(node u) except +
		count degreeOut(node u) except +
		node addNode() except +
		void removeNode(node u) except +
		void addEdge(node u, node v, edgeweight w) except +
		void removeEdge(node u, node v) except +
		bool hasEdge(node u, node v) except +
		edgeweight weight(node u, node v) except +
		vector[node] nodes() except +
		vector[pair[node, node]] edges() except +
		vector[node] neighbors(node u) except +
		bool isWeighted() except +
		bool isDirected() except +
		string toString() except +
		string getName() except +
		void setName(string name) except +
		edgeweight totalEdgeWeight() except +
		node randomNode() except +
		node randomNeighbor(node) except +
		pair[node, node] randomEdge() except +


cdef class Graph:
	""" An undirected graph (with optional weights) and parallel iterator methods. 
		
		Graph(n=0, weighted=False, directed=False)

		Create a graph of `n` nodes. The graph has assignable edge weights if `weighted` is set to True.
	 	If `weighted` is set to False each edge has edge weight 1.0 and any other weight assignment will 
	 	be ignored.

	    Parameters
	    ----------
	    n : count, optional
	    	Number of nodes.
	    weighted : bool, optional
	    	If set to True, the graph can have edge weights other than 1.0.
	    directed : bool, optional
	    	If set to True, the graph will be directed.
	"""
	cdef _Graph* _this

	def __cinit__(self, n=0, weighted=False, directed=False):
		self._this = new _Graph(n, weighted, directed)

	# # any _thisect which appears as a return type needs to implement setThis
	# cdef setThis(self, _Graph other):
	# 	#del self._this
	# 	self._this = other
	# 	return self

	cdef setThis(self, _Graph* other):
		self._this = other
		return self

	# this is necessary so that the C++ object gets properly garbage collected
	def __dealloc__(self):
		del self._this

	def numberOfNodes(self):
		""" 
		Get the number of nodes in the graph.
	 	
	 	Returns
	 	-------
	 	count
	 		The number of nodes.
		"""
		return self._this.numberOfNodes()

	def numberOfEdges(self):
		""" 
		Get the number of edges in the graph.
	 	
	 	Returns
	 	-------
	 	count
	 		The number of edges.
		"""
		return self._this.numberOfEdges()

	def degree(self, u):
		""" 
		Get the number of neighbors of `v`.
	 	
		Parameters
		----------
		v : node
			Node.

		Returns
		-------
		count
			The number of neighbors.
		"""
		return self._this.degree(u)

	def degreeIn(self, u):
		return self._this.degreeIn(u)

	def degreeOut(self, u):
		return self._this.degreeOut(u)

	def addNode(self):
		""" Add a new node to the graph and return it.

		Returns
		-------
		node
			The new node.
	 	"""
		return self._this.addNode()

	def removeNode(self, u):
		""" Remove the isolated node `u` from the graph.
	 	
	 	Parameters
	 	----------
	 	u : node
	 		Node.

	 	Notes
	 	-----
	 	Although it would be convenient to remove all incident edges at the same time, this causes complications for 
	 	dynamic applications. Therefore, removeNode is an atomic event. All incident edges need to be removed first 
	 	and an exception is thrown otherwise.
		"""
		self._this.removeNode(u)

	def addEdge(self, u, v, w=1.0):
		""" Insert an undirected edge between the nodes `u` and `v`. If the graph is weighted you can optionally
	 	set a weight for this edge. The default weight is 1.0.

	 	Parameters
	 	----------
	 	u : node
	 		Endpoint of edge.
 		v : node
 			Endpoint of edge.
		w : edgeweight, optional
			Edge weight.
		"""
		self._this.addEdge(u, v, w)

	def removeEdge(self, u, v):
		""" Removes the undirected edge {`u`,`v`}.

		Parameters
		----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.
		"""
		self._this.removeEdge(u, v)

	def hasEdge(self, u, v):
		""" Checks if undirected edge {`u`,`v`} exists in the graph.

		Parameters
		----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.

		Returns
		-------
		bool
			True if the edge exists, False otherwise.
		"""
		return self._this.hasEdge(u, v)

	def weight(self, u, v):
		""" Get edge weight of edge {`u` , `v`}. Returns 0 if edge does not exist.

		Parameters
		----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.

		Returns
		-------
		edgeweight
			Edge weight of edge {`u` , `v`} or 0 if edge does not exist.
		"""
		return self._this.weight(u, v)

	def nodes(self):
		""" Get list of all nodes.
	 	
	 	Returns
	 	-------
	 	list
	 		List of all nodes.
		"""
		return self._this.nodes()

	def edges(self):
		""" Get list of edges as node pairs.
	 
	 	Returns
	 	-------
	 	list
	 		List of edges as node pairs.
		"""
		return self._this.edges()

	def neighbors(self, u):
		""" Get list of neighbors of `u`.
	 	
	 	Parameters
	 	----------
	 	u : node
	 		Node.

	 	Returns 
	 	-------
	 	list
	 		List of neighbors of `u.
		"""
		return self._this.neighbors(u)

	def isWeighted(self):
		"""
		Returns
		-------
		bool
			True if this graph supports edge weights other than 1.0.	 
		"""
		return self._this.isWeighted()

	def isDirected(self):
		return self._this.isDirected()

	def toString(self):
		""" Get a string representation of the graph.

		Returns
		-------
		string
			A string representation of the graph.
		"""
		return self._this.toString()

	def getName(self):
		""" Get the name of the graph.

		Returns
		-------
		string
			The name of the graph.
		"""
		return pystring(self._this.getName())

	def setName(self, name):
		""" Set name of graph to `name`.

		Parameters
		----------
		name : string
			The name.
		"""
		self._this.setName(stdstring(name))

	def totalEdgeWeight(self):
		""" Get the sum of all edge weights.

		Returns
		-------
		edgeweight
			The sum of all edge weights.
		"""
		return self._this.totalEdgeWeight()

	def randomNode(self):
		""" Get a random node of the graph.

		Returns
		-------
		node
			A random node.
		"""
		return self._this.randomNode()

	def randomNeighbor(self, u):
		""" Get a random neighbor of `v` and `none` if degree is zero.

		Parameters
		----------
		v : node
			Node.

		Returns
		-------
		node
			A random neighbor of `v.
		"""
		return self._this.randomNeighbor(u)

	def randomEdge(self):
		""" Get a random edge of the graph.

		Returns
		-------
		pair
			Random random edge.

		Notes
		-----
		Fast, but not uniformly random.
		"""
		return self._this.randomEdge()


# TODO: expose all methods

cdef extern from "../cpp/graph/BFS.h":
	cdef cppclass _BFS "NetworKit::BFS":
		_BFS(_Graph G, node source) except +
		void run() except +
		vector[edgeweight] getDistances() except +
		vector[node] getPath(node t) except +

cdef class BFS:
	""" Simple breadth-first search on a Graph from a given source

	BFS(G, source)

	Create BFS for `G` and source node `source`.

	Parameters
	----------
	G : Graph
		The graph.
	source : node
		The source node of the breadth-first search.		

	"""
	cdef _BFS* _this

	def __cinit__(self, Graph G, source):		
		self._this = new _BFS(dereference(G._this), source)

	def run(self):
		"""	
		Breadth-first search from source.
		
		Returns
		-------
		vector
			Vector of unweighted distances from source node, i.e. the
	 		length (number of edges) of the shortest path from source to any other node.
		"""
		self._this.run()

	def getDistances(self):
		""" 
		Returns a vector of weighted distances from the source node, i.e. the
 	 	length of the shortest path from the source node to any other node.

 	 	Returns
 	 	-------
 	 	vector
 	 		The weighted distances from the source node to any other node in the graph.
		"""
		return self._this.getDistances()

	def getPath(self, t):
		""" Returns a shortest path from source to `t` and an empty path if source and `t` are not connected.

		Parameters
		----------
		t : node
			Target node.

		Returns
		-------
		vector
			A shortest path from source to `t or an empty path.
		"""
		return self._this.getPath(t)


cdef extern from "../cpp/graph/Dijkstra.h":
	cdef cppclass _Dijkstra "NetworKit::Dijkstra":
		_Dijkstra(_Graph G, node source) except +
		void run() except +
		vector[edgeweight] getDistances() except +
		vector[node] getPath(node t) except +

cdef class Dijkstra:
	""" Dijkstra's SSSP algorithm.
	Returns list of weighted distances from node source, i.e. the length of the shortest path from source to 
	any other node.

    Dijkstra(G, source)

    Creates Dijkstra for `G` and source node `source`.

    Parameters
	----------
	G : Graph
		The graph.
	source : node
		The source node.
    """
	cdef _Dijkstra* _this

	def __cinit__(self, Graph G, source):		
		self._this = new _Dijkstra(dereference(G._this), source)

	def run(self):
		""" Performs the Dijkstra SSSP algorithm on the graph given in the constructor. """
		self._this.run()

	def getDistances(self):
		""" Returns a vector of weighted distances from the source node, i.e. the
 	 	length of the shortest path from the source node to any other node.

 	 	Returns
 	 	-------
 	 	vector
 	 		The weighted distances from the source node to any other node in the graph.
		"""
		return self._this.getDistances()

	def getPath(self, t):
		""" Returns a shortest path from source to `t` and an empty path if source and `t` are not connected.

		Parameters
		----------
		t : node
			Target node.

		Returns
		-------
		vector
			A shortest path from source to `t or an empty path.
		"""
		return self._this.getPath(t)


cdef extern from "../cpp/graph/Subgraph.h":
	cdef cppclass _Subgraph "NetworKit::Subgraph":
		_Subgraph() except +
		_Graph* _fromNodes(_Graph G, unordered_set[node] nodes)

cdef class Subgraph:
	""" Methods for creating subgraphs """
	cdef _Subgraph _this

	def fromNodes(self, Graph G, nodes): #unordered_set[node]
		""" Create a subgraph induced by the set `nodes`.
	 	
	 	Parameters
	 	----------
	 	G : Graph
	 		The graph.
 		nodes : list
 			A subset of nodes of `G` which induce the subgraph.

		Returns
		-------
		Graph
			The subgraph induced by `nodes`.

		Notes
		-----
		The returned graph G' is isomorphic (structurally identical) to the subgraph in G,
	 	but node indices are not preserved.
		"""
		cdef unordered_set[node] nnodes
		for node in nodes:
			nnodes.insert(node);
		return Graph().setThis(self._this._fromNodes(dereference(G._this), nnodes))

cdef extern from "../cpp/independentset/Luby.h":
	cdef cppclass _Luby "NetworKit::Luby":
		_Luby() except +
		vector[bool] run(_Graph G)
		string toString()


# FIXME: check correctness
cdef class Luby:
	""" Luby's parallel maximal independent set algorithm"""
	cdef _Luby _this

	def run(self, Graph G not None):
		""" Returns a boolean vector of length n where vec[v] is True iff v is in the independent sets.

		Parameters
		----------
		G : Graph
			The graph.

		Returns
		-------
		vector
			A boolean vector of length n.
		"""
		return self._this.run(dereference(G._this))

	def toString(self):
		""" Get string representation of the algorithm.

		Returns
		-------
		string
			The string representation of the algorithm.
		"""
		return self._this.toString().decode("utf-8")


# Module: generators

# cdef extern from "../cpp/graph/GraphGenerator.h":
# 	cdef cppclass _GraphGenerator "NetworKit::GraphGenerator":
# 		_GraphGenerator() except +
# 		_Graph makeRandomGraph(count n, double p)


# cdef class GraphGenerator:
# 	""" Provides several functions for graph generation"""
# 	cdef _GraphGenerator _this

# 	def __cinit__(self):
# 		self._this = _GraphGenerator()


# 	def makeRandomGraph(self, n, p):
# 		cdef _Graph _G = self._this.makeRandomGraph(n, p)
# 		return Graph(0).setThis(_G)

cdef extern from "../cpp/generators/BarabasiAlbertGenerator.h":
	cdef cppclass _BarabasiAlbertGenerator "NetworKit::BarabasiAlbertGenerator":
		_BarabasiAlbertGenerator() except +
		_BarabasiAlbertGenerator(count k, count nMax, count n0) except +
		#_Graph* _generate()
		_Graph* _generate() except +

cdef class BarabasiAlbertGenerator:
	""" Generates a scale-free graph using the Barabasi-Albert preferential attachment model. """
	cdef _BarabasiAlbertGenerator _this

	def __cinit__(self, k, nMax, n0):
		""" TODO
		"""
		self._this = _BarabasiAlbertGenerator(k, nMax, n0)

	def generate(self):
		""" TODO
		"""
		return Graph().setThis(self._this._generate());


cdef extern from "../cpp/generators/PubWebGenerator.h":
	cdef cppclass _PubWebGenerator "NetworKit::PubWebGenerator":
		_PubWebGenerator(count numNodes, count numberOfDenseAreas, float neighborhoodRadius, count maxNumberOfNeighbors) except +
		_Graph* _generate() except +

cdef class PubWebGenerator:
	""" Generates a static graph that resembles an assumed geometric distribution of nodes in 
	a P2P network. 

	The basic structure is to distribute points randomly in the unit torus 
	and to connect vertices close to each other (at most @a neighRad distance and none of 
	them already has @a maxNeigh neighbors). The distribution is chosen to get some areas with 
	high density and others with low density. There are @a numDenseAreas dense areas, which can 
	overlap. Each area is circular, has a certain position and radius and number of points. 
	These values are strored in @a denseAreaXYR and @a numPerArea, respectively.

	Used and described in more detail in J. Gehweiler, H. Meyerhenke: A Distributed
	Diffusive Heuristic for Clustering a Virtual P2P Supercomputer. In Proc. 7th High-Performance
	Grid Computing Workshop (HPGC'10), in conjunction with 24th IEEE Internatl. Parallel and
	Distributed Processing Symposium (IPDPS'10), IEEE, 2010.

	PubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	Parameters
	----------
	numNodes : count
		Up to a few thousand (possibly more if visualization is not desired and quadratic 
		time complexity has been resolved)
	numberOfDenseAreas : count
		Depending on number of nodes, e.g. [8, 50]
	neighborhoodRadius : float
		The higher, the better the connectivity [0.1, 0.35]
	maxNumberOfNeighbors : count
		Maximum degree, a higher value corresponds to better connectivity [4, 40]
	"""
	cdef _PubWebGenerator* _this

	def __cinit__(self, numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors):
		""" TODO
		"""
		self._this = new _PubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	def generate(self):
		""" TODO
		"""
		return Graph(0).setThis(self._this._generate())


cdef extern from "../cpp/generators/ErdosRenyiGenerator.h":
	cdef cppclass _ErdosRenyiGenerator "NetworKit::ErdosRenyiGenerator":
		_ErdosRenyiGenerator(count nNodes, double prob) except +
		_Graph* _generate() except +

cdef class ErdosRenyiGenerator:
	""" Creates random graphs in the G(n,p) model.
	The generation follows Vladimir Batagelj and Ulrik Brandes: "Efficient
	generation of large random networks", Phys Rev E 71, 036113 (2005).	 

	ErdosRenyiGenerator(count, double)

	Creates G(nNodes, prob) graphs.

	Parameters
	----------
	nNodes : count
		Number of nodes n in the graph.
	prob : double
		Probability of existence for each edge p.		
	"""

	cdef _ErdosRenyiGenerator* _this

	def __cinit__(self, nNodes, prob):		
		self._this = new _ErdosRenyiGenerator(nNodes, prob)

	def generate(self):
		return Graph(0).setThis(self._this._generate())


cdef extern from "../cpp/generators/DorogovtsevMendesGenerator.h":
	cdef cppclass _DorogovtsevMendesGenerator "NetworKit::DorogovtsevMendesGenerator":
		_DorogovtsevMendesGenerator(count nNodes) except +
		_Graph* _generate() except +

cdef class DorogovtsevMendesGenerator:
	"""
	TODO:

	Parameters
	----------
	nNodes : count
		Number of nodes in the target graph.
	"""

	cdef _DorogovtsevMendesGenerator* _this

	def __cinit__(self, nNodes):
		self._this = new _DorogovtsevMendesGenerator(nNodes)

	def generate(self):
		return Graph(0).setThis(self._this._generate())


cdef extern from "../cpp/generators/ClusteredRandomGraphGenerator.h":
	cdef cppclass _ClusteredRandomGraphGenerator "NetworKit::ClusteredRandomGraphGenerator":
		_ClusteredRandomGraphGenerator(count, count, double, double) except +
		_Graph* _generate() except +

cdef class ClusteredRandomGraphGenerator:
	""" The ClusteredRandomGraphGenerator class is used to create a clustered random graph.
 		
	The number of nodes and the number of edges are adjustable as well as the probabilities
	for intra-cluster and inter-cluster edges.

	ClusteredRandomGraphGenerator(count, count, pin, pout)

	Creates a clustered random graph.

	Parameters
	----------
	n : count
		number of nodes
	k : count
		number of clusters
	pin : double
		intra-cluster edge probability
	pout : double
		inter-cluster edge probability
	"""

	cdef _ClusteredRandomGraphGenerator* _this

	def __cinit__(self, n, k, pin, pout):		
		self._this = new _ClusteredRandomGraphGenerator(n, k, pin, pout)

	def generate(self):
		""" Generates a clustered random graph with the properties given in the constructor.

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph(0).setThis(self._this._generate())


cdef extern from "../cpp/generators/ChungLuGenerator.h":
	cdef cppclass _ChungLuGenerator "NetworKit::ChungLuGenerator":
		# TODO: revert to count when cython issue fixed
		_ChungLuGenerator(vector[unsigned int] degreeSequence) except +
		_Graph* _generate() except +

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
		""" Generates graph with expected degree sequence seq. 

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph(0).setThis(self._this._generate())


cdef extern from "../cpp/generators/HavelHakimiGenerator.h":
	cdef cppclass _HavelHakimiGenerator "NetworKit::HavelHakimiGenerator":
		# TODO: revert to count when cython issue fixed
		_HavelHakimiGenerator(vector[unsigned int] degreeSequence, bool skipTest) except +
		_Graph* _generate() except +
		bool isRealizable() except +
		bool getRealizable() except +

cdef class HavelHakimiGenerator:
	""" Havel-Hakimi algorithm for generating a graph according to a given degree sequence.
 		
 		The sequence, if it is realizable, is reconstructed exactly. The resulting graph usually 
 		has a high clustering coefficient. Construction runs in linear time O(m). However, the test 
 		if a sequence is realizable is quadratic in the sequence length.

 		HavelHakimiGenerator(sequence, skipTest=True)

 		Parameters
		----------
		sequence : vector
			Degree sequence to realize. Must be non-increasing.
		skipTest : bool, optional
			If True, the test if the sequence is realizable is skipped.
	        Default value is False. Set ONLY to True if you are certain that the
	        sequence is realizable		
	"""

	cdef _HavelHakimiGenerator* _this


	def __cinit__(self, degreeSequence, skipTest=True):		
		self._this = new _HavelHakimiGenerator(degreeSequence, skipTest)

	def isRealizable(self):
		return self._this.isRealizable()

	def getRealizable(self):
		return self._this.getRealizable();

	def generate(self):
		""" Generates degree sequence seq (if it is realizable).

		Returns
		-------
		Graph
			Empty graph if graph is not realizable, otherwise graph with degree sequence seq.
		"""
		return Graph(0).setThis(self._this._generate())


cdef extern from "../cpp/generators/RmatGenerator.h":
	cdef cppclass _RmatGenerator "NetworKit::RmatGenerator":
		_RmatGenerator(count scale, count edgeFactor, double a, double b, double c, double d) except +
		_Graph* _generate() except +

cdef class RmatGenerator:
	"""
	Generates static R-MAT graphs. R-MAT (recursive matrix) graphs are
	random graphs with n=2^scale nodes and m=nedgeFactor edges.
	More details at http://www.graph500.org or in the original paper:
	Deepayan Chakrabarti, Yiping Zhan, Christos Faloutsos:
	R-MAT: A Recursive Model for Graph Mining. SDM 2004: 442-446.	

	RmatGenerator(scale, edgeFactor, a, b, c, d)

	Parameters
	----------
	scale : count
		Number of nodes = 2^scale
	edgeFactor : count
		Number of edges = number of nodes * edgeFactor
	a : double
		Probability for quadrant upper left
	b : double
		Probability for quadrant upper right
	c : double
		Probability for quadrant lower left
	d : double
		Probability for quadrant lower right
	
	"""

	cdef _RmatGenerator* _this

	def __cinit__(self, count scale, count edgeFactor, double a, double b, double c, double d):		
		self._this = new _RmatGenerator(scale, edgeFactor, a, b, c, d)

	def generate(self):
		""" Graph to be generated according to parameters specified in constructor.

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph(0).setThis(self._this._generate())


# Module: graphio

cdef extern from "../cpp/io/METISGraphReader.h":
	cdef cppclass _METISGraphReader "NetworKit::METISGraphReader":
		_METISGraphReader() except +
		# _Graph read(string path) except +
		_Graph* _read(string path) except +
		_Graph* readToHeap(string path) except +

cdef class METISGraphReader:
	""" Reads the METIS adjacency file format [1]. If the Fast reader fails,
		use readGraph(path, graphio.formats.metis) as an alternative.
		[1]: http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
	"""
	cdef _METISGraphReader _this

	def read(self, path):
		pathbytes = path.encode("utf-8") # string needs to be converted to bytes, which are coerced to std::string
		return Graph(0).setThis(self._this._read(pathbytes))


cdef extern from "../cpp/io/EdgeListReader.h":
	cdef cppclass _EdgeListReader "NetworKit::EdgeListReader":
		_EdgeListReader() except +
		_EdgeListReader(char separator, node firstNode, string commentPrefix, bool continuous)
		_Graph read(string path) except +
		_Graph* _read(string path) except +
		unordered_map[node,node] getNodeMap() except +
		

cdef class EdgeListReader:
	""" Reads the METIS adjacency file format [1]. If the Fast reader fails,
		use readGraph(path, graphio.formats.metis) as an alternative.
		[1]: http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
	"""
	cdef _EdgeListReader _this

	def __cinit__(self, separator, firstNode, commentPrefix="#", continuous=True):
		self._this = _EdgeListReader(stdstring(separator)[0], firstNode, stdstring(commentPrefix), continuous)

	def read(self, path):
		pathbytes = path.encode("utf-8") # string needs to be converted to bytes, which are coerced to std::string
		return Graph(0).setThis(self._this._read(pathbytes))
	
	def getNodeMap(self):
		cdef unordered_map[node,node] cResult = self._this.getNodeMap()
		result = []
		for elem in cResult:
			result.append((elem.first,elem.second))
		return result

cdef extern from "../cpp/io/METISGraphWriter.h":
	cdef cppclass _METISGraphWriter "NetworKit::METISGraphWriter":
		_METISGraphWriter() except +
		void write(_Graph G, string path) except +


cdef class METISGraphWriter:
	""" Writes graphs in the METIS format"""
	cdef _METISGraphWriter _this

	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		self._this.write(dereference(G._this), stdstring(path))


cdef extern from "../cpp/io/DotGraphWriter.h":
	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter":
		_DotGraphWriter() except +
		void write(_Graph G, string path) except +


cdef class DotGraphWriter:
	""" Writes graphs in the .dot/GraphViz format"""
	cdef _DotGraphWriter _this

	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		self._this.write(dereference(G._this), stdstring(path))


#cdef extern from "../cpp/io/VNAGraphWriter.h":
#	cdef cppclass _VNAGraphWriter "NetworKit::VNAGraphWriter":
#		_VNAGraphWriter() except +
#		void write(_Graph G, string path) except +


#cdef class VNAGraphWriter:
#	""" Writes graphs in the VNA format. The VNA format is commonly used by Netdraw, and is very similar to Pajek format.
#	It defines nodes and edges (ties), and supports attributes. Each section of the file is separated by an asterisk. """
#	cdef _VNAGraphWriter _this

#	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
#		self._this.write(dereference(G._this), stdstring(path))


cdef extern from "../cpp/io/GMLGraphWriter.h":
	cdef cppclass _GMLGraphWriter "NetworKit::GMLGraphWriter":
		_GMLGraphWriter() except +
		void write(_Graph G, string path) except +


cdef class GMLGraphWriter:
	""" Writes a graph and its coordinates as a GML file.[1]
		[1] http://svn.bigcat.unimaas.nl/pvplugins/GML/trunk/docs/gml-technical-report.pdf """
	cdef _GMLGraphWriter _this

	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		self._this.write(dereference(G._this), stdstring(path))


cdef extern from "../cpp/io/EdgeListWriter.h":
	cdef cppclass _EdgeListWriter "NetworKit::EdgeListWriter":
		_EdgeListWriter() except +
		_EdgeListWriter(char separator, node firstNode) except +
		void write(_Graph G, string path) except +

cdef class EdgeListWriter:
	""" Reads and writes graphs in various edge list formats. The constructor takes a
		seperator char and the ID of the first node as paraneters."""

	cdef _EdgeListWriter _this

	def __cinit__(self, separator, firstNode):
		cdef char sep = stdstring(separator)[0]
		self._this = _EdgeListWriter(sep, firstNode)

	def write(self, Graph G not None, path):
		self._this.write(dereference(G._this), stdstring(path))



cdef extern from "../cpp/io/LineFileReader.h":
	cdef cppclass _LineFileReader "NetworKit::LineFileReader":
		_LineFileReader() except +
		vector[string] read(string path)


cdef class LineFileReader:
	""" Reads a file and puts each line in a list of strings """
	cdef _LineFileReader _this

	def read(self, path):
		return self._this.read(stdstring(path))


cdef extern from "../cpp/io/SNAPGraphWriter.h":
	cdef cppclass _SNAPGraphWriter "NetworKit::SNAPGraphWriter":
		_SNAPGraphWriter() except +
		void write(_Graph G, string path)

cdef class SNAPGraphWriter:
	""" Writes graphs in a format suitable for the Georgia Tech SNAP software [1]
		[1]: http://snap-graph.sourceforge.net/
	"""
	cdef _SNAPGraphWriter _this

	def write(self, Graph G, path):
		self._this.write(dereference(G._this), stdstring(path))


cdef extern from "../cpp/io/SNAPGraphReader.h":
	cdef cppclass _SNAPGraphReader "NetworKit::SNAPGraphReader":
		_SNAPGraphReader() except +
		_Graph read(string path) except +
		_Graph* _read(string path) except +
		unordered_map[node,node] getNodeIdMap() except +

cdef class SNAPGraphReader:
	""" Reads a graph from the SNAP graph data collection [1] (currently experimental)
		[1]: http://snap.stanford.edu/data/index.html
	"""
	cdef _SNAPGraphReader _this

	def read(self, path):
		return Graph().setThis(self._this._read(stdstring(path)))

	def getNodeIdMap(self):
		cdef unordered_map[node,node] cResult = self._this.getNodeIdMap()
		result = []
		for elem in cResult:
			result.append((elem.first,elem.second))
		return result


cdef extern from "../cpp/io/PartitionReader.h":
	cdef cppclass _PartitionReader "NetworKit::PartitionReader":
		_PartitionReader() except +
		_Partition read(string path)


cdef class PartitionReader:
	""" Reads a partition from a file.
		File format: line i contains subset id of element i.
	 """
	cdef _PartitionReader _this

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))


cdef extern from "../cpp/io/PartitionWriter.h":
	cdef cppclass _PartitionWriter "NetworKit::PartitionWriter":
		_PartitionWriter() except +
		void write(_Partition, string path)


cdef class PartitionWriter:
	""" Writes a partition to a file.
		File format: line i contains subset id of element i.
	 """
	cdef _PartitionWriter _this

	def write(self, Partition zeta, path):
		self._this.write(zeta._this, stdstring(path))


cdef extern from "../cpp/io/EdgeListPartitionReader.h":
	cdef cppclass _EdgeListPartitionReader "NetworKit::EdgeListPartitionReader":
		_EdgeListPartitionReader() except +
		_EdgeListPartitionReader(node firstNode) except +
		_Partition read(string path)


cdef class EdgeListPartitionReader:
	""" Reads a partition from an edge list type of file
	 """
	cdef _EdgeListPartitionReader _this

	def __cinit__(self, firstNode=1):
		self._this = _EdgeListPartitionReader(firstNode)

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))

cdef extern from "../cpp/io/SNAPEdgeListPartitionReader.h":
	cdef cppclass _SNAPEdgeListPartitionReader "NetworKit::SNAPEdgeListPartitionReader":
		_SNAPEdgeListPartitionReader() except +
		_Cover read(string path, unordered_map[node,node] nodeMap,_Graph G) except +
#		_Partition readWithInfo(string path, count nNodes) except +

cdef class SNAPEdgeListPartitionReader:
	""" Reads a partition from a SNAP 'community with ground truth' file
	 """
	cdef _SNAPEdgeListPartitionReader _this

	def read(self,path, nodeMap, Graph G):
		cdef unordered_map[node,node] cNodeMap
		for (key,val) in nodeMap:
			cNodeMap[key] = val
		return Cover().setThis(self._this.read(stdstring(path), cNodeMap, dereference(G._this)))

#	def readWithInfo(self,path,nNodes):
#		return Partition().setThis(self._this.readWithInfo(stdstring(path),nNodes))

#not existing yet, maybe in the future?
#cdef extern from "../cpp/io/EdgeListPartitionWriter.h":
#	cdef cppclass _EdgeListPartitionWriter "NetworKit::EdgeListPartitionWriter":
#		_EdgeListPartitionWriter() except +
#		void write(_Partition, string path)


#cdef class EdgeListPartitionWriter:
#	""" Writes a partition to a edge list type of file.
#		File format: a line contains the element id and the subsed id of the element.
#	 """
#	cdef _EdgeListPartitionWriter _this

#	def Write(self, Partition zeta, path):
#		self._this.write(zeta._this, stdstring(path))

cdef extern from "../cpp/io/CoverReader.h":
	cdef cppclass _CoverReader "NetworKit::CoverReader":
		_CoverReader() except +
		_Cover read(string path,_Graph G) except +

cdef class CoverReader:
	""" Reads a cover from a file
		File format: each line contains the space-separated node ids of a community
	 """
	cdef _CoverReader _this

	def read(self, path, Graph G):
		return Cover().setThis(self._this.read(stdstring(path), dereference(G._this)))

cdef extern from "../cpp/io/CoverWriter.h":
	cdef cppclass _CoverWriter "NetworKit::CoverWriter":
		_CoverWriter() except +
		void write(_Cover, string path)


cdef class CoverWriter:
	""" Writes a partition to a file.
		File format: each line contains the space-separated node ids of a community
	 """
	cdef _CoverWriter _this

	def write(self, Cover zeta, path):
		self._this.write(zeta._this, stdstring(path))

cdef extern from "../cpp/io/EdgeListCoverReader.h":
	cdef cppclass _EdgeListCoverReader "NetworKit::EdgeListCoverReader":
		_EdgeListCoverReader() except +
		_EdgeListCoverReader(node firstNode) except +
		_Cover read(string path, _Graph G) except +


cdef class EdgeListCoverReader:
	""" Reads a cover from an edge list type of file
		File format: each line starts with a node id and continues with a list of the communities the node belongs to
	 """
	cdef _EdgeListCoverReader _this

	def __cinit__(self, firstNode=1):
		self._this = _EdgeListCoverReader(firstNode)

	def read(self, path, Graph G):
		return Cover().setThis(self._this.read(stdstring(path), dereference(G._this)))

# Parameters

cdef extern from "../cpp/base/Parameters.h":
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


# Module: structures
#
cdef extern from "../cpp/structures/Partition.h":
	cdef cppclass _Partition "NetworKit::Partition":
		_Partition() except +
		_Partition(index) except +
		index subsetOf(index e) except +
		index extend() except +
		void remove(index e) except +
		void addToSubset(index s, index e) except +
		void moveToSubset(index s, index e) except +
		void toSingleton(index e) except +
		void allToSingletons() except +
		void mergeSubsets(index s, index t) except +
		void setUpperBound(index upper) except +
		index upperBound() except +
		index lowerBound() except +
		void compact() except +
		bool contains(index e) except +
		bool inSameSubset(index e1, index e2) except +
		vector[count] subsetSizes() except +
		map[index, count] subsetSizeMap() except +
		set[index] getMembers(const index s) except +
		count numberOfElements() except +
		count numberOfSubsets() except +
		vector[index] getVector() except +
		void setName(string name) except +
		string getName() except +
		set[index] getSubsetIds() except +


cdef class Partition:
	""" Implements a partition of a set, i.e. a subdivision of the 
 		set into disjoint subsets.

 		Partition(z=0)

 		Create a new partition data structure for `z` elements.

		Parameters
		----------
		z : index, optional
			Maximum index of an element. Default is 0.		
	"""
	cdef _Partition _this

	def __cinit__(self, z=0):		
		self._this = _Partition(z)

	def __len__(self):
		"""
		Returns
		-------
		count
			Number of elements in the partition.
		"""
		return self._this.numberOfElements()

	def __getitem__(self, e):
		""" Get the set (id) in which the element `e` is contained.
	 
	 	Parameters
	 	----------
	 	e : index
	 		Index of element.

	 	Returns
	 	-------
	 	index
	 		The index of the set in which `e` is contained.
		"""
		return self._this.subsetOf(e)

	cdef setThis(self, _Partition other):
		self._this = other
		return self

	def subsetOf(self, e):
		""" Get the set (id) in which the element `e` is contained.
	 
	 	Parameters
	 	----------
	 	e : index
	 		Index of element.

	 	Returns
	 	-------
	 	index
	 		The index of the set in which `e` is contained.
		"""
		return self._this.subsetOf(e)

	def extend(self):
		""" Extend the data structure and create a slot	for one more element.

		Initializes the entry to `none` and returns the index of the entry.

		Returns
		-------
		index
			The index of the new element.
		"""
		self._this.extend()

	def addToSubset(self, s, e):
		""" Add a (previously unassigned) element `e` to the set `s`.

		Parameters
		----------
		s : index
			The index of the subset.
		e : index
			The element to add.
		"""
		self._this.addToSubset(s, e)

	def moveToSubset(self, index s, index e):
		"""  Move the (previously assigned) element `e` to the set `s.

		Parameters
		----------
		s : index
			The index of the subset.
		e : index
			The element to move.
		"""
		self._this.moveToSubset(s, e)

	def toSingleton(self, index e):
		""" Creates a singleton set containing the element `e`.

		Parameters
		----------
		e : index
			The index of the element.
		"""
		self._this.toSingleton(e)

	def allToSingletons(self):
		""" Assigns every element to a singleton set. Set id is equal to element id. """
		self._this.allToSingletons()

	def mergeSubsets(self, index s, index t):
		""" Assigns the elements from both sets to a new set and returns the id of it.

		Parameters
		----------
		s : index
			Set to merge.
		t : index
			Set to merge.

		Returns
		-------
		index
			Id of newly created set.
		"""
		self._this.mergeSubsets(s, t)


	def setUpperBound(self, index upper):
		""" Sets an upper bound for the subset ids that **can** be assigned.

		Parameters
		----------
		upper : index
			Highest assigned subset id + 1
		"""
		self._this.setUpperBound(upper)

	def upperBound(self):
		""" Return an upper bound for the subset ids that have been assigned.
	 	(This is the maximum id + 1.)

	 	Returns
	 	-------
	 	index
	 		The upper bound.
		"""
		return self._this.upperBound()

	def lowerBound(self):
		""" Get a lower bound for the subset ids that have been assigned.

		Returns
		-------
		index
			The lower bound.
		"""
		return self._this.lowerBound()

	def compact(self):
		""" Change subset IDs to be consecutive, starting at 0. """
		self._this.compact()

	def contains(self, index e):
		""" Check if partition assigns a valid subset to the element `e`.

		Parameters
		----------
		e : index
			The element.

		Returns
		-------
		bool
			True if the assigned subset is valid, False otherwise.
		"""
		return self._this.contains(e)

	def inSameSubset(self, index e1, index e2):
		""" Check if two elements `e1` and `e2` belong to the same subset.

		Parameters
		----------
		e1 : index
			An Element.
		e2 : index
			An Element.

		Returns
		-------
		bool
			True if `e1` and `e2` belong to same subset, False otherwise.
		"""
		return self._this.inSameSubset(e1, e2)

	def subsetSizes(self):
		""" Get a list of subset sizes. Indices do not necessarily correspond to subset ids.
	 	
	 	Returns
	 	-------
	 	vector
	 		A vector of subset sizes.
		"""
		return self._this.subsetSizes()

	def subsetSizeMap(self):
		""" Get a map from subset id to size of the subset.

		Returns
		-------
		dict
			A map from subset id to size of the subset.
		"""
		return self._this.subsetSizeMap()

	def getMembers(self, s):
		""" Get the members of the subset `s`.

		Parameters
		----------
		s : index
			The subset.

		Returns
		-------
		set
			A set containing the members of `s.
		"""
		return self._this.getMembers(s)

	def numberOfElements(self):
		"""
		Returns
		-------
		count
			Number of elements in the partition.
		"""
		return self._this.numberOfElements()

	def numberOfSubsets(self):
		""" Get the current number of sets in this partition.

		Returns
		-------
		count
			The current number of sets.
		"""
		return self._this.numberOfSubsets()

	def getVector(self):
		""" Get the actual vector representing the partition data structure.

		Returns
		-------
		vector
			Vector containing information about partitions.
		"""
		return self._this.getVector()

	def setName(self, string name):
		"""  Set a human-readable identifier `name` for the instance.

		Parameters
		----------
		name : string
			The name.
		"""
		self._this.setName(name)

	def getName(self):
		""" Get the human-readable identifier.

		Returns
		-------
		string
			The name of this partition.
		"""
		return self._this.getName()

	def getSubsetIds(self):
		""" Get the ids of nonempty subsets.

		Returns
		-------
		set
			A set of ids of nonempty subsets.
		"""
		return self._this.getSubsetIds()


cdef extern from "../cpp/structures/Cover.h":
	cdef cppclass _Cover "NetworKit::Cover":
		_Cover() except +
		set[index] subsetsOf(index e) except +
#		index extend() except +
		void remove(index e) except +
		void addToSubset(index s, index e) except +
		void moveToSubset(index s, index e) except +
		void toSingleton(index e) except +
		void allToSingletons() except +
		void mergeSubsets(index s, index t) except +
#		void setUpperBound(index upper) except +
		index upperBound() except +
		index lowerBound() except +
#		void compact() except +
		bool contains(index e) except +
		bool inSameSubset(index e1, index e2) except +
		vector[count] subsetSizes() except +
		map[index, count] subsetSizeMap() except +
		set[index] getMembers(const index s) except +
		count numberOfElements() except +
		count numberOfSubsets() except +
#		vector[index] getVector() except +
#		void setName(string name) except +
#		string getName() except +
#		set[index] getSubsetIds() except +


cdef class Cover:
	""" Implements a cover of a set, i.e. an assignment of its elements to possibly overlapping subsets. """
	cdef _Cover _this

	cdef setThis(self, _Cover other):
		self._this = other
		return self

	def subsetsOf(self, e):
		""" Get the ids of subsets in which the element `e` is contained.

		Parameters
		----------
		e : index
			An element

		Returns
		-------
		set
			A set of subset ids in which `e` 	is contained.
		"""
		return self._this.subsetsOf(e)

#	def extend(self):
#		self._this.extend()

	def addToSubset(self, s, e):
		""" Add the (previously unassigned) element `e` to the set `s`.

		Parameters
		----------
		s : index
			A subset
		e : index
			An element			
		"""
		self._this.addToSubset(s, e)

	def moveToSubset(self, index s, index e):
		""" Move the element `e` to subset `s`, i.e. remove it from all other subsets and place it in the subset.

		Parameters
		----------
		s : index
			A subset
		e : index
			An element
		"""
		self._this.moveToSubset(s, e)

	def toSingleton(self, index e):
		""" Creates a singleton set containing the element `e` and returns the index of the new set.

		Parameters
		----------
		e : index
			An element

		Returns
		-------
		index
			The index of the new set.
		"""
		self._this.toSingleton(e)

	def allToSingletons(self):
		""" Assigns every element to a singleton set. Set id is equal to element id. """
		self._this.allToSingletons()

	def mergeSubsets(self, index s, index t):
		""" Assigns the elements from both sets to a new set.

		Parameters
		----------
		s : index
			A subset
		t : index
			A subset
		"""
		self._this.mergeSubsets(s, t)

#	def setUpperBound(self, index upper):
#		self._this.setUpperBound(upper)

	def upperBound(self):
		""" Get an upper bound for the subset ids that have been assigned.
	   	(This is the maximum id + 1.)

	   	Returns
	   	-------
	   	index
	   		An upper bound.
		"""
		return self._this.upperBound()

	def lowerBound(self):
		""" Get a lower bound for the subset ids that have been assigned.

		Returns
		-------
		index
			A lower bound.
		"""
		return self._this.lowerBound()

#	def compact(self):
#		self._this.compact()

	def contains(self, index e):
		"""  Check if cover assigns a valid subset to the element `e`.

		Parameters
		----------
		e : index
			An element.

		Returns
		-------
		bool
			True, if `e` is assigned to a valid subset, False otherwise.

		"""
		return self._this.contains(e)

	def inSameSubset(self, index e1, index e2):
		"""  Check if two elements `e1` and `e2` belong to the same subset.
	 
	 	Parameters
	 	----------
	 	e1 : index
			An element.
		e2 : index
			An element.

		Returns
		-------
		bool
			True, if `e1` and `e2` belong to the same subset, False otherwise.
		"""
		return self._this.inSameSubset(e1, e2)

	def subsetSizes(self):
		""" Get a list of subset sizes. 

		Returns
		-------
		list
			A list of subset sizes.

		Notes
		-----
		Indices do not necessarily correspond to subset ids.
		"""
		return self._this.subsetSizes()

	def subsetSizeMap(self):
		""" Get a map from subset id to size of the subset.
	 
	 	Returns
	 	-------
	 	dict
	 		A map from subset id to size of the subset.
		"""
		return self._this.subsetSizeMap()

	def getMembers(self, s):
		""" Get the members of a specific subset `s`.

		Returns
		-------
		set
			The set of members of subset `s`.
		"""
		return self._this.getMembers(s)

	def numberOfElements(self):
		""" Get the current number of elements in this cover.

		Returns
		-------
		count
			The current number of elements.
		"""
		return self._this.numberOfElements()

	def numberOfSubsets(self):
		"""  Get the current number of sets in this cover.

		Returns
		-------
		count
			The number of sets in this cover.
		"""
		return self._this.numberOfSubsets()

#	def getVector(self):
#		return self._this.getVector()

#	def setName(self, string name):
#		self._this.setName(name)

#	def getName(self):
#		return self._this.getName()

#	def getSubsetIds(self):
#		return self._this.getSubsetIds()


# Module: community

cdef extern from "../cpp/community/Coverage.h":
	cdef cppclass _Coverage "NetworKit::Coverage":
		_Coverage() except +
		double getQuality(_Partition _zeta, _Graph _G) except +

cdef class Coverage:
	""" Coverage is the fraction of intra-community edges """
	cdef _Coverage _this

	def getQuality(self, Partition zeta, Graph G):
		return self._this.getQuality(zeta._this, dereference(G._this))


cdef extern from "../cpp/community/Modularity.h":
	cdef cppclass _Modularity "NetworKit::Modularity":
		_Modularity() except +
		double getQuality(_Partition _zeta, _Graph _G) except +


cdef class Modularity:
	"""	Modularity is a quality index for community detection. 
	It assigns a quality value in [-0.5, 1.0] to a partition of a graph which is higher for more modular networks and 
	partitions which better capture the modular structure. See also http://en.wikipedia.org/wiki/Modularity_(networks).

 	Notes
	-----
	Modularity is defined as:

	.. math:: mod(\zeta) := \\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)} - \\frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }

	"""
	cdef _Modularity _this

	def getQuality(self, Partition zeta, Graph G):
		return self._this.getQuality(zeta._this, dereference(G._this))


cdef class CommunityDetector:
	""" Abstract base class for static community detection algorithms """
	pass

cdef extern from "../cpp/community/PLP.h":
	cdef cppclass _PLP "NetworKit::PLP":
		_PLP() except +
		_PLP(count updateThreshold) except +
		_Partition run(_Graph _G) except +
		_Partition runFromGiven(_Graph _G, _Partition _part) except +
		count numberOfIterations() except +
		string toString() except +


cdef class PLP(CommunityDetector):
	""" Parallel label propagation for community detection:
	Moderate solution quality, very short time to solution.

	Notes
	-----
	As described in Ovelgoenne et al: An Ensemble Learning Strategy for Graph Clustering
 	Raghavan et al. proposed a label propagation algorithm for graph clustering.
 	This algorithm initializes every vertex of a graph with a unique label. Then, in iterative
 	sweeps over the set of vertices the vertex labels are updated. A vertex gets the label
 	that the maximum number of its neighbors have. The procedure is stopped when every vertex
 	has the label that at least half of its neighbors have.
	"""
	cdef _PLP _this

	def __cinit__(self, updateThreshold=None):
		if updateThreshold is None:
			self._this = _PLP()
		else:
			self._this = _PLP(updateThreshold)


	def run(self, Graph G not None):
		""" Run the label propagation clustering algorithm.
		
		Parameters
		----------
		G : Graph
			input graph

	 	Returns
	 	-------
	 	Partition
	 		The created clustering.
		"""
		return Partition().setThis(self._this.run(dereference(G._this)))

	def runFromGiven(self, Graph G not None, Partition part not None):
		""" Run the label propagation clustering algorithm starting
		from the Partition part.
		
		Parameters
		----------
		G : Graph
			input graph
			
		part : Partition
			input partition

	 	Returns
	 	-------
	 	Partition
	 		The created clustering.
		"""
		return Partition().setThis(self._this.run(dereference(G._this)))

	def numberOfIterations(self):
		""" Get number of iterations in last run.

		Returns
		-------
		count
			The number of iterations.
		"""
		return self._this.numberOfIterations()

	def toString(self):
		""" Get string representation.

		Returns
		-------
		string
			String representation of algorithm and parameters.
		"""
		return self._this.toString().decode("utf-8")


cdef extern from "../cpp/community/LPDegreeOrdered.h":
	cdef cppclass _LPDegreeOrdered "NetworKit::LPDegreeOrdered":
		_LPDegreeOrdered() except +
		_Partition run(_Graph _G)
		count numberOfIterations()

cdef class LPDegreeOrdered(CommunityDetector):
	""" Label propagation-based community detection algorithm which processes nodes in increasing order of node degree.	"""
	cdef _LPDegreeOrdered _this

	def run(self, Graph G not None):
		return Partition().setThis(self._this.run(dereference(G._this)))

	def numberOfIterations(self):
		""" Get number of iterations in last run.

		Returns
		-------
		count
			Number of iterations.
		"""
		return self._this.numberOfIterations()


cdef extern from "../cpp/community/PLM.h":
	cdef cppclass _PLM "NetworKit::PLM":
		_PLM() except +
		_PLM(bool refine, double gamma, string par, count maxIter) except +
		string toString() except +
		_Partition run(_Graph G) except +


cdef class PLM(CommunityDetector):
	""" MultiLevel Parallel LocalMover - the Louvain method, optionally extended to
		a full multi-level algorithm with refinement

		PLM(refine=True, gamma=1.0, par="balanced", maxIter=32)

		Parameters
		----------
		refine : bool, optional
			Add a second move phase to refine the communities.
		gamma : double
			Multi-resolution modularity parameter:
			1.0 -> standard modularity
	 		0.0 -> one community
	 		2m 	-> singleton communities
		par : string
			parallelization strategy
		maxIter : count
			maximum number of iterations for move phase
	"""

	cdef _PLM _this

	def __cinit__(self, refine=True, gamma=1.0, par="balanced", maxIter=32):
		self._this = _PLM(refine, gamma, stdstring(par), maxIter)

	def toString(self):
		""" Get string representation.

		Returns
		-------
		string
			String representation of this algorithm.
		"""
		return self._this.toString().decode("utf-8")

	def run(self, Graph G not None):
		""" Detect communities in the given graph `G`

		Parameters
		----------
		G : Graph
			The graph.

		Returns
		-------
		Partition
			A partition containing the found communities.
		"""
		return Partition().setThis(self._this.run(dereference(G._this)))


cdef extern from "../cpp/community/CNM.h":
	cdef cppclass _CNM "NetworKit::CNM":
		string toString() except +
		_Partition run(_Graph G) except +


cdef class CNM(CommunityDetector):
	"""
	Community detection algorithm due to Clauset, Newman and Moore.
 	Probably not the fastest possible implementation, but it already uses a priority queue
 	and local updates.
 	"""

	cdef _CNM* _this

	def __cinit__(self):
		self._this = new _CNM()

	def toString(self):
		""" Get string representation.

		Returns
		-------
		string
			A string representation of this algorithm.
		"""
		return self._this.toString().decode("utf-8")

	def run(self, Graph G not None):
		""" Detect communities in the given graph `graph`.

		Parameters
		----------
		graph : Graph
			The graph.
		
		Returns
		-------
		Partition
			A partition containing the found communities.
		"""
		return Partition().setThis(self._this.run(dereference(G._this)))


cdef class DissimilarityMeasure:
	""" Abstract base class for partition/community dissimilarity measures """
	pass


cdef extern from "../cpp/community/NodeStructuralRandMeasure.h":
	cdef cppclass _NodeStructuralRandMeasure "NetworKit::NodeStructuralRandMeasure":
		_NodeStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second)

cdef class NodeStructuralRandMeasure(DissimilarityMeasure):
	""" The node-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering all pairs of nodes.
	"""
	cdef _NodeStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		return self._this.getDissimilarity(dereference(G._this), first._this, second._this)


cdef extern from "../cpp/community/GraphStructuralRandMeasure.h":
	cdef cppclass _GraphStructuralRandMeasure "NetworKit::GraphStructuralRandMeasure":
		_GraphStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second)

cdef class GraphStructuralRandMeasure(DissimilarityMeasure):
	""" The graph-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering connected pairs of nodes.
	"""
	cdef _GraphStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		return self._this.getDissimilarity(dereference(G._this), first._this, second._this)


cdef extern from "../cpp/community/JaccardMeasure.h":
	cdef cppclass _JaccardMeasure "NetworKit::JaccardMeasure":
		_JaccardMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second)

cdef class JaccardMeasure(DissimilarityMeasure):
	""" TODO:
	"""
	cdef _JaccardMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		return self._this.getDissimilarity(dereference(G._this), first._this, second._this)

cdef extern from "../cpp/community/NMIDistance.h":
	cdef cppclass _NMIDistance "NetworKit::NMIDistance":
		_NMIDistance() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second)

cdef class NMIDistance(DissimilarityMeasure):
	""" The NMI distance assigns a similarity value in [0,1] to two partitions
		of a graph.
	"""
	cdef _NMIDistance _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		return self._this.getDissimilarity(dereference(G._this), first._this, second._this)

cdef extern from "../cpp/community/EPP.h":
	cdef cppclass _EPP "NetworKit::EPP":
		_Partition run(_Graph G)
		string toString()

cdef class EPP(CommunityDetector):
	""" EPP - Ensemble Preprocessing community detection algorithm.
	Combines multiple base algorithms and a final algorithm. A consensus of the 
	solutions of the base algorithms is formed and the graph is coarsened accordingly. 
	Then the final algorithm operates on the coarse graph and determines a solution 
	for the input graph.
	"""
	cdef _EPP _this

	def run(self, Graph G):
		"""  Run the ensemble clusterer on `G` and return the result in a Partition.

		Parameters
		----------
		G : Graph
			The graph.

		Returns
		-------
		Partition:
			A Partition of the clustering.
		"""
		return Partition().setThis(self._this.run(dereference(G._this)))

	def toString(self):
		""" String representation of EPP class.

		Returns
		-------
		string
			String representation.
		"""
		return self._this.toString()

	cdef setThis(self, _EPP other):
		self._this = other
		return self


cdef extern from "../cpp/community/EPPFactory.h":
	cdef cppclass _EPPFactory "NetworKit::EPPFactory":
		_EPP make(count ensembleSize, string baseAlgorithm, string finalAlgorithm)

cdef class EPPFactory:
	""" This class makes instaces of the EPP community detection algorithm """
	cdef _EPPFactory _this

	def make(self, ensembleSize, baseAlgorithm="PLP", finalAlgorithm="PLM"):
		return EPP().setThis(self._this.make(ensembleSize, stdstring(baseAlgorithm), stdstring(finalAlgorithm)))

cdef extern from "../cpp/community/CommunityGraph.h":
	cdef cppclass _CommunityGraph "NetworKit::CommunityGraph":
		void run(_Graph G, _Partition zeta) except +
		_Graph* _getGraph() except +
		map[index, node] getCommunityToNodeMap() except +
		map[node, index] getNodeToCommunityMap() except +

cdef class CommunityGraph:
	""" The CommunityGraph class represents a Graph coarsened according to communities. Each node in the CommunityGraph 
 	represents a community. Edge weights are the weights of inter-community cuts.
	"""
	cdef _CommunityGraph _this

	def run(self, Graph G, Partition zeta):
		""" Creates a coarsened graph of `G` according to communities in `zeta`. Edge weights are the weights of 
		inter-community cuts.

		Parameters
		----------
		G : Graph
			The graph.
		zeta : Partition
			A community clustering of `G`.
		"""
		self._this.run(dereference(G._this), zeta._this)

	def getGraph(self):
		""" Returns the coarsened Graph.

		Returns
		-------
		Graph
			The coarsened graph.
		"""
		return Graph().setThis(self._this._getGraph())

	def getCommunityToNodeMap(self):
		""" Maps community id to node id in the community graph.

		Returns
		-------
		dict
			Map containing community id to node id mappings.
		"""
		return self._this.getCommunityToNodeMap()

	def getNodeToCommunityMap(self):
		""" Maps node id in the community graph to community id.

		Returns
		-------
		dict
			Map containing node id to community id mappins.
		"""
		return self._this.getNodeToCommunityMap()

# Module: properties

# this is an example for using static methods
cdef extern from "../cpp/properties/GraphProperties.h" namespace "NetworKit::GraphProperties":
	# static methods live in the class namespace, so declare them here
	pair[count, count] minMaxDegree(_Graph _G) except +
	double averageDegree(_Graph _G) except +
	vector[count] degreeDistribution(_Graph _G) except +
	vector[double] localClusteringCoefficients(_Graph _G) except +
	double averageLocalClusteringCoefficient(_Graph _G) except +
	vector[double] localClusteringCoefficientPerDegree(_Graph _G) except +
	double degreeAssortativity(_Graph G, bool) except +

	cdef cppclass _GraphProperties "NetworKit::GraphProperties":
		pass

cdef class GraphProperties:
	""" Collects various functions for basic graph properties """

	@staticmethod
	def minMaxDegree(Graph G not None):
		return minMaxDegree(dereference(G._this))

	@staticmethod
	def averageDegree(Graph G not None):
		return averageDegree(dereference(G._this))

	@staticmethod
	def degreeDistribution(Graph G not None):
		return degreeDistribution(dereference(G._this))

	@staticmethod
	def averageLocalClusteringCoefficient(Graph G not None):
		""" The average local clustering coefficient for the graph `G`.

		Parameters
		----------
		G : Graph
			The graph.

		Notes
		-----

		.. math:: \\frac{1}{n} \cdot \sum_{v \in V} c_v

		"""
		return averageLocalClusteringCoefficient(dereference(G._this))

	@staticmethod
	def degreeAssortativity(Graph G, bool useWeights):	
		""" Get degree assortativity of the graph `G`.

		Parameters
		----------
		G : Graph
			The graph
		useWeights : bool
			If True, the weights are considered for calculation.

		Returns
		-------
		double
			Degree assortativity of the graph `G`.

		Notes
		-----
		Degree assortativity based on description in Newman: Networks. An Introduction. Chapter 8.7.
		"""
		return degreeAssortativity(dereference(G._this), useWeights)




cdef extern from "../cpp/properties/ConnectedComponents.h":
	cdef cppclass _ConnectedComponents "NetworKit::ConnectedComponents":
		_ConnectedComponents(_Graph G) except +
		void run() except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		_Partition getPartition() except +
		map[index, count] getComponentSizes() except +


cdef class ConnectedComponents:
	""" Determines the connected components and associated values for an undirected graph.

	ConnectedComponents(G)

	Create ConnectedComponents for Graph `G`.

	Parameters
	----------
	G : Graph
		The graph.
	"""
	cdef _ConnectedComponents* _this

	def __cinit__(self,  Graph G):
		self._this = new _ConnectedComponents(dereference(G._this))

	def run(self):
		""" This method determines the connected components for the graph given in the constructor. """
		self._this.run()

	def getPartition(self):
		""" Get a Partition that represents the components.

		Returns
		-------
		Partition
			A partition representing the found components.
		"""
		return Partition().setThis(self._this.getPartition())

	def numberOfComponents(self):
		""" Get the number of connected components.

		Returns
		-------
		count:
			The number of connected components.
		"""
		return self._this.numberOfComponents()

	def componentOfNode(self, v):
		"""  Get the the component in which node `v` is situated.

		v : node
			The node whose component is asked for.
		"""
		return self._this.componentOfNode(v)

	def getComponentSizes(self):
		return self._this.getComponentSizes()


cdef extern from "../cpp/properties/ParallelConnectedComponents.h":
	cdef cppclass _ParallelConnectedComponents "NetworKit::ParallelConnectedComponents":
		_ParallelConnectedComponents(_Graph G, bool coarsening) except +
		void run() except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		_Partition getPartition() except +


cdef class ParallelConnectedComponents:
	""" Determines the connected components and associated values for
		an undirected graph.
	"""
	cdef _ParallelConnectedComponents* _this

	def __cinit__(self,  Graph G, coarsening=True	):
		self._this = new _ParallelConnectedComponents(dereference(G._this), coarsening)

	def run(self):
		self._this.run()

	def getPartition(self):
		return Partition().setThis(self._this.getPartition())

	def numberOfComponents(self):
		return self._this.numberOfComponents()

	def componentOfNode(self, v):
		return self._this.componentOfNode(v)


cdef extern from "../cpp/properties/StronglyConnectedComponents.h":
	cdef cppclass _StronglyConnectedComponents "NetworKit::StronglyConnectedComponents":
		_StronglyConnectedComponents(_Graph G) except +
		void run() except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		_Partition getPartition() except +


cdef class StronglyConnectedComponents:
	""" Determines the connected components and associated values for
		a directed graph.
	"""
	cdef _StronglyConnectedComponents* _this

	def __cinit__(self,  Graph G):
		self._this = new _StronglyConnectedComponents(dereference(G._this))

	def run(self):
		self._this.run()

	def getPartition(self):
		return Partition().setThis(self._this.getPartition())

	def numberOfComponents(self):
		return self._this.numberOfComponents()

	def componentOfNode(self, v):
		return self._this.componentOfNode(v)



cdef extern from "../cpp/properties/ClusteringCoefficient.h" namespace "NetworKit::ClusteringCoefficient":
		vector[double] exactLocal(_Graph G) except +
		double avgLocal(_Graph G) except +
		double approxAvgLocal(_Graph G, count trials) except +
		double exactGlobal(_Graph G) except +
		double approxGlobal(_Graph G, count trials) except +

cdef class ClusteringCoefficient:
	@staticmethod
	def exactLocal(Graph G):
		return exactLocal(dereference(G._this))

	@staticmethod	
	def avgLocal(Graph G):
		"""  This calculates the average local clustering coefficient of graph `G`.

		Parameters
		----------
		G : Graph
			The graph.

		Notes
		-----

		.. math:: c(G) := \\frac{1}{n} \sum_{u \in V} c(u)

		where

		.. math:: c(u) := \\frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}

		"""
		return avgLocal(dereference(G._this))

	@staticmethod
	def approxAvgLocal(Graph G, trials):
		return approxAvgLocal(dereference(G._this), trials)

	@staticmethod
	def exactGlobal(Graph G):
		""" This calculates the global clustering coefficient. """
		return exactGlobal(dereference(G._this))

	@staticmethod
	def approxGlobal(Graph G, trials):
		return approxGlobal(dereference(G._this), trials)



cdef extern from "../cpp/properties/Diameter.h" namespace "NetworKit::Diameter":
	pair[count, count] estimatedDiameterRange(_Graph G, double error) except +
	count exactDiameter(_Graph G) except +
	edgeweight estimatedVertexDiameter(_Graph G, count) except +

cdef class Diameter:
	"""
	TODO: docstring
	"""

	@staticmethod
	def estimatedDiameterRange(Graph G, error=0.1):
		""" Estimates a range for the diameter of @a G. Based on the algorithm suggested in 
		C. Magnien, M. Latapy, M. Habib: Fast Computation of Empirically Tight Bounds for 
		the Diameter of Massive Graphs. Journal of Experimental Algorithmics, Volume 13, Feb 2009.

		Returns
		-------
		pair
			Pair of lower and upper bound for diameter.
		"""
		return estimatedDiameterRange(dereference(G._this), error)

	@staticmethod
	def exactDiameter(Graph G):
		""" Get the exact diameter of the graph `G`.

		Parameters
		----------
		G : Graph
			The graph.

		Returns
		-------
		edgeweight
			Exact diameter of the graph `G`.
		"""
		return exactDiameter(dereference(G._this))

	@staticmethod
	def estimatedVertexDiameter(Graph G, samples):
		""" Get a 2-approximation of the node diameter (unweighted diameter) of `G`.

		Parameters
		----------
		G : Graph
			The graph.
		samples : count
			One sample is enough if the graph is connected. If there 
			are multiple connected components, then the number of samples 
			must be chosen so that the probability of sampling the component 
			with the largest diameter ist high. 

		Returns
		-------
		edgeweight
			A 2-approximation of the vertex diameter (unweighted diameter) of `G`.
		"""
		return estimatedVertexDiameter(dereference(G._this), samples)

cdef extern from "../cpp/properties/Eccentricity.h" namespace "NetworKit::Eccentricity":
	pair[node, edgeweight] getValue(_Graph G, node v) except +

cdef class Eccentricity:
	"""
	TODO: docstring
	"""

	@staticmethod
	def getValue(Graph G, v):
		return getValue(dereference(G._this), v)


cdef extern from "../cpp/properties/CoreDecomposition.h":
	cdef cppclass _CoreDecomposition "NetworKit::CoreDecomposition":
		_CoreDecomposition(_Graph)
		void run() except +
		vector[index] coreNumbers() except +
		index coreNumber(node) except +
		vector[set[node]] cores() except +
		vector[set[node]] shells() except +
		index maxCoreNumber() except +

cdef class CoreDecomposition:
	""" Computes k-core decomposition of a graph.

	CoreDecomposition(G)

	Create CoreDecomposition class for graph `G`.

	Parameters
	----------
	G : Graph
		The graph.
	"""

	cdef _CoreDecomposition* _this

	def __cinit__(self, Graph G):
		self._this = new _CoreDecomposition(dereference(G._this))

	def run(self):
		""" Perform k-core decomposition of graph passed in constructor. """
		self._this.run()

	def coreNumbers(self):
		""" Get vector of core numbers, indexed by node.

		Returns
		-------
		vector
			Vector of core numbers, indexed by node.
		"""
		return self._this.coreNumbers()

	def coreNumber(self, v):
		""" Get core number of node `v`.

		Parameters
		----------
		v : node
			A node.

		Returns
		-------
		node
			Core number of node `v.
		"""
		return self._this.coreNumber(v)

	def maxCoreNumber(self):
		""" Get maximum core number.

		Returns
		-------
		index
			The maximum core number.
		"""
		return self._this.maxCoreNumber()

	def cores(self):
		""" Get the k-cores as sets of nodes, indexed by k.

		Returns
		-------
		vector
			The k-cores as sets of nodes, indexed by k.
		"""
		return self._this.cores()

	def shells(self):
		""" Get the k-shells as sets of nodes, indexed by k.

		Returns
		-------
		vector
			The k-shells as sets of nodes, indexed by k.
		"""
		return self._this.shells()


# Module: centrality


# TODO: how to properly wrap class hierarchies and reuse code?

# cdef extern from "../cpp/centrality/Centrality.h":
# 	cdef cppclass _Centrality "NetworKit::Centrality":
# 		_centrality(_Graph, bool) except +
# 		void run() except +
# 		vector[double] scores() except +
# 		vector[pair[node, double]] ranking() except +
# 		double score(node) except +


# cdef class Centrality:
# 	""" Abstract base class for centrality measures"""

# 	def __cinit__(self, _Centrality* _this):
# 		self._this = _this

# 	def run(self):
# 		self._this.run()

# 	def scores(self):
# 		return self._this.scores()

# 	def score(self, v):
# 		return self._this.score(v)

# 	def ranking(self):
# 		return self._this.ranking()


cdef extern from "../cpp/centrality/Betweenness.h":
	cdef cppclass _Betweenness "NetworKit::Betweenness":
		_Betweenness(_Graph, bool) except +
		void run() except +
		void run(bool runUnweightedInParallel) except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +

cdef class Betweenness:
	"""
		Betweenness(G, normalized=False)

		Constructs the Betweenness class for the given Graph `G`. If the betweenness scores should be normalized,
  		then set `normalized` to True.
	 
	 	Parameters
	 	----------
	 	G : Graph
	 		The graph.
	 	normalized : bool, optional
	 		Set this parameter to True if scores should be normalized in the interval [0,1].
	"""
	cdef _Betweenness* _this

	def __cinit__(self, Graph G, normalized=False):
		self._this = new _Betweenness(dereference(G._this), normalized)

	def run(self, parallel=False):
		"""  Compute betweenness scores sequential or parallel depending on `runUnweightedInParallel`.

		Parameters
		----------
		runUnweightedInParallel : bool
			If set to True the computation is done in parallel.
		"""
		self._this.run(parallel)

	def scores(self):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""
		return self._this.scores()

	def score(self, v):
		""" Get the betweenness score of node `v` calculated by run().

		Parameters
		----------
		v : node
			A node.

		Returns
		-------
		double
			The betweenness score of node `v.
		"""
		return self._this.score(v)

	def ranking(self):
		""" Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score 
		calculated by run().

		Returns
		-------
		vector
			A vector of pairs.
		"""
		return self._this.ranking()


cdef extern from "../cpp/centrality/ApproxBetweenness.h":
	cdef cppclass _ApproxBetweenness "NetworKit::ApproxBetweenness":
		_ApproxBetweenness(_Graph, double, double, count) except +
		void run() except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +
		count numberOfSamples() except +

cdef class ApproxBetweenness:
	""" Approximation of betweenness centrality according to algorithm described in 
 	Matteo Riondato and Evgenios M. Kornaropoulos: Fast Approximation of Betweenness Centrality through Sampling

 	ApproxBetweenness(G, epsiolon=0.01, delta=0.1)

 	The algorithm approximates the betweenness of all vertices so that the scores are
	within an additive error epsilon with probability at least (1- delta).
	The values are normalized by default.

	Parameters
	----------
	G : Graph
		the graph
	epsilon : double, optional
		maximum additive error
	delta : double, optional
		probability that the values are within the error guarantee
	"""
	cdef _ApproxBetweenness* _this

	def __cinit__(self, Graph G, epsilon=0.01, delta=0.1, diameterSamples=0):
		self._this = new _ApproxBetweenness(dereference(G._this), epsilon, delta, diameterSamples)

	def run(self):
		self._this.run()

	def scores(self):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""
		return self._this.scores()

	def score(self, v):
		""" Get the betweenness score of node `v` calculated by run().

		Parameters
		----------
		v : node
			A node.

		Returns
		-------
		double
			The betweenness score of node `v.
		"""
		return self._this.score(v)

	def ranking(self):
		""" Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score 
		calculated by run().

		Returns
		-------
		vector
			A vector of pairs.
		"""
		return self._this.ranking()

	def numberOfSamples(self):
		return self._this.numberOfSamples()


cdef extern from "../cpp/centrality/ApproxBetweenness2.h":
	cdef cppclass _ApproxBetweenness2 "NetworKit::ApproxBetweenness2":
		_ApproxBetweenness2(_Graph, count, bool) except +
		void run() except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +

cdef class ApproxBetweenness2:
	""" Approximation of betweenness centrality according to algorithm described in 
	Sanders, Geisberger, Schultes: Better Approximation of Betweenness Centrality

	ApproxBetweenness2(G, nSamples, normalized=False)

	The algorithm approximates the betweenness of all nodes, using weighting 
	of the contributions to avoid biased estimation.

	Parameters
	----------
	G : Graph
		input graph
	nSamples : count
		user defined number of samples
	normalized : bool, optional
		normalize centrality values in interval [0,1]
	"""
	cdef _ApproxBetweenness2* _this

	def __cinit__(self, Graph G, nSamples, normalized=False):
		self._this = new _ApproxBetweenness2(dereference(G._this), nSamples, normalized)

	def run(self):
		self._this.run()

	def scores(self):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""
		return self._this.scores()

	def score(self, v):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""
		return self._this.score(v)

	def ranking(self):
		""" Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score 
		calculated by run().

		Returns
		-------
		vector
			A vector of pairs.
		"""
		return self._this.ranking()


cdef extern from "../cpp/centrality/PageRank.h":
	cdef cppclass _PageRank "NetworKit::PageRank":
		_PageRank(_Graph, double damp, double tol) except +
		void run() except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +

cdef class PageRank:
	"""	Compute PageRank as node centrality measure.

	PageRank(G, damp, tol=1e-9)

	Parameters
	----------
	G : Graph
		Graph to be processed.
	damp : double
		Damping factor of the PageRank algorithm.
	tol : double, optional
		Error tolerance for PageRank iteration.
	"""
	cdef _PageRank* _this

	def __cinit__(self, Graph G, double damp, double tol=1e-9):
		self._this = new _PageRank(dereference(G._this), damp, tol)

	def run(self):
		self._this.run()

	def scores(self):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""		
		return self._this.scores()

	def score(self, v):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""
		return self._this.score(v)

	def ranking(self):
		""" Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score 
		calculated by run().

		Returns
		-------
		vector
			A vector of pairs.
		"""
		return self._this.ranking()


cdef extern from "../cpp/centrality/EigenvectorCentrality.h":
	cdef cppclass _EigenvectorCentrality "NetworKit::EigenvectorCentrality":
		_EigenvectorCentrality(_Graph, double tol) except +
		void run() except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +

cdef class EigenvectorCentrality:
	"""	Computes the leading eigenvector of the graph's adjacency matrix (normalized in 2-norm). 
	Interpreted as eigenvector centrality score.

	EigenvectorCentrality(G, tol=1e-9)

	Constructs the EigenvectorCentrality class for the given Graph `G`. `tol` defines the tolerance for convergence.

	Parameters
	----------
	G : Graph
		The graph.
	tol : double, optional
		The tolerance for convergence.
	"""
	cdef _EigenvectorCentrality* _this

	def __cinit__(self, Graph G, double tol=1e-9):
		self._this = new _EigenvectorCentrality(dereference(G._this), tol)

	def run(self):
		self._this.run()

	def scores(self):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""		
		return self._this.scores()

	def score(self, v):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""
		return self._this.score(v)

	def ranking(self):
		""" Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score 
		calculated by run().

		Returns
		-------
		vector
			A vector of pairs.
		"""
		return self._this.ranking()


cdef extern from "../cpp/centrality/DegreeCentrality.h":
	cdef cppclass _DegreeCentrality "NetworKit::DegreeCentrality":
		_DegreeCentrality(_Graph, bool normalized) except +
		void run() except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +

cdef class DegreeCentrality:
	""" Node centrality index which ranks nodes by their degree.
 	Optional normalization by maximum degree.

 	DegreeCentrality(G, normalized=False)

 	Constructs the DegreeCentrality class for the given Graph `G`. If the betweenness scores should be normalized, 
 	then set `normalized` to True.

 	Parameters
 	----------
 	G : Graph
 		The graph.
 	normalized : bool, optional
 		Normalize centrality values in the interval [0,1].
	"""
	cdef _DegreeCentrality* _this

	def __cinit__(self, Graph G, bool normalized=False):
		self._this = new _DegreeCentrality(dereference(G._this), normalized)

	def run(self):
		self._this.run()

	def scores(self):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""		
		return self._this.scores()

	def score(self, v):
		""" Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""
		return self._this.score(v)

	def ranking(self):
		""" Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score 
		calculated by run().

		Returns
		-------
		vector
			A vector of pairs.
		"""
		return self._this.ranking()


# Module: dynamic

cdef extern from "../cpp/dynamics/GraphEvent.h":
	enum _Type "NetworKit::GraphEvent::Type":
		NODE_ADDITION,
		NODE_REMOVAL,
		EDGE_ADDITION,
		EDGE_REMOVAL,
		EDGE_WEIGHT_UPDATE,
		TIME_STEP

cdef extern from "../cpp/dynamics/GraphEvent.h":
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
			self._this.v = v

	property w:
		def __get__(self):
			return self._this.w
		def __set__(self, w):
			self._this.w = w

	def __cinit__(self, type, u, v, w):
		self._this = _GraphEvent(type, u, v, w)

	def toString(self):
		return self._this.toString().decode("utf-8")

	def __repr__(self):
		return self.toString()


cdef extern from "../cpp/dynamics/DGSStreamParser.h":
	cdef cppclass _DGSStreamParser "NetworKit::DGSStreamParser":
		_DGSStreamParser(string path, bool mapped, node baseIndex) except +
		vector[_GraphEvent] getStream() except +

cdef class DGSStreamParser:
	cdef _DGSStreamParser* _this

	def __cinit__(self, path, mapped=True, baseIndex=0):
		self._this = new _DGSStreamParser(stdstring(path), mapped, baseIndex)

	def getStream(self):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.getStream()]


cdef extern from "../cpp/dynamics/DGSWriter.h":
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


# cdef extern from "../cpp/dcd2/DynamicCommunityDetection.h":
# 	cdef cppclass _DynamicCommunityDetection "NetworKit::DynamicCommunityDetection":
# 		_DynamicCommunityDetection(string inputPath, string algoName, string updateStrategy, count interval, count restart, vector[string] recordSettings) except +
# 		void run() except +
# 		vector[double] getTimeline(string key) except +
# 		vector[pair[count, count]] getGraphSizeTimeline() except +
# 		vector[pair[_Graph, _Partition]] getResultTimeline() except +

# cdef class DynamicCommunityDetection:
# 	cdef _DynamicCommunityDetection* _this

# 	def __cinit__(self, inputPath, algoName, updateStrategy, interval, restart, recordSettings):
# 		self._this = new _DynamicCommunityDetection(stdstring(inputPath), stdstring(algoName), stdstring(updateStrategy), interval, restart, [stdstring(key) for key in recordSettings])

# 	def run(self):
# 		self._this.run()

# 	def getTimeline(self, key):
# 		return self._this.getTimeline(stdstring(key))

# 	def getGraphSizeTimeline(self):
# 		return self._this.getGraphSizeTimeline()

# 	def getResultTimeline(self):
# 		timeline = []
# 		for pair in self._this.getResultTimeline():
# 			_G = pair.first
# 			_zeta = pair.second
# 			timeline.append((Graph().setThis(_G), Partition().setThis(_zeta)))
# 		return timeline



cdef extern from "../cpp/generators/DynamicPathGenerator.h":
	cdef cppclass _DynamicPathGenerator "NetworKit::DynamicPathGenerator":
		_DynamicPathGenerator() except +
		vector[_GraphEvent] generate(count nSteps) except +


cdef class DynamicPathGenerator:
	""" Example dynamic graph generator: Generates a dynamically growing path. """
	cdef _DynamicPathGenerator* _this

	def __cinit__(self):
		self._this = new _DynamicPathGenerator()

	def generate(self, nSteps):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]


cdef extern from "../cpp/generators/DynamicDorogovtsevMendesGenerator.h":
	cdef cppclass _DynamicDorogovtsevMendesGenerator "NetworKit::DynamicDorogovtsevMendesGenerator":
		_DynamicDorogovtsevMendesGenerator() except +
		vector[_GraphEvent] generate(count nSteps) except +


cdef class DynamicDorogovtsevMendesGenerator:
	cdef _DynamicDorogovtsevMendesGenerator* _this

	def __cinit__(self):
		self._this = new _DynamicDorogovtsevMendesGenerator()

	def generate(self, nSteps):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]



cdef extern from "../cpp/generators/DynamicPubWebGenerator.h":
	cdef cppclass _DynamicPubWebGenerator "NetworKit::DynamicPubWebGenerator":
		_DynamicPubWebGenerator(count numNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors) except +
		vector[_GraphEvent] generate(count nSteps) except +
		_Graph* _getGraph() except +


cdef class DynamicPubWebGenerator:
	cdef _DynamicPubWebGenerator* _this

	def __cinit__(self, numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors):
		self._this = new _DynamicPubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	def generate(self, nSteps):
		""" Generate event stream.

		Parameters
		----------
		nSteps : count
			Number of time steps in the event stream.
		"""
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]

	def getGraph(self):
		return Graph().setThis(self._this._getGraph())


# cdef extern from "../cpp/generators/ForestFireGenerator.h":
# 	cdef cppclass _ForestFireGenerator "NetworKit::ForestFireGenerator":
# 		_ForestFireGenerator(double p) except +
# 		vector[_GraphEvent] generate(count nSteps) except +
# 		_Graph getGraph() except +


# cdef class ForestFireGenerator:
# 	cdef _ForestFireGenerator* _this

# 	def __cinit__(self, p):
# 		self._this = new _ForestFireGenerator(p)

# 	def generate(self, nSteps):
# 		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]



cdef extern from "../cpp/dynamics/GraphUpdater.h":
	cdef cppclass _GraphUpdater "NetworKit::GraphUpdater":
		_GraphUpdater(_Graph G) except +
		void update(vector[_GraphEvent] stream) except +
		vector[pair[count, count]] getSizeTimeline() except +

cdef class GraphUpdater:
	cdef _GraphUpdater* _this

	def __cinit__(self, Graph G):
		self._this = new _GraphUpdater(dereference(G._this))

	def update(self, stream):
		cdef vector[_GraphEvent] _stream
		for ev in stream:
			_stream.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.update(_stream)
