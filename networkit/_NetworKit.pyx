
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
from networkit.unordered_set cimport unordered_set
from networkit.unordered_map cimport unordered_map

# NetworKit typedefs
ctypedef uint64_t count
ctypedef uint64_t index
ctypedef uint64_t edgeid
ctypedef index node
ctypedef index cluster
ctypedef double edgeweight

cdef extern from "cpp/Globals.h" namespace "NetworKit":
	index _none "NetworKit::none"

none = _none

cdef extern from "<algorithm>" namespace "std":
	void swap[T](T &a,  T &b)
	_Graph move( _Graph t ) nogil # specialized declaration as general declaration disables template argument deduction and doesn't work
	_Partition move( _Partition t) nogil
	_Cover move(_Cover t) nogil
	pair[_Graph, vector[node]] move(pair[_Graph, vector[node]]) nogil

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

# Cython helper functions

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes

def pystring(stdstring):
	""" convert a std::string (= python byte string) to a normal Python string"""
	return stdstring.decode("utf-8")


# Function definitions

cdef extern from "cpp/auxiliary/Log.h" namespace "Aux":
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

cdef extern from "cpp/auxiliary/Parallelism.h" namespace "Aux":
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


cdef extern from "cpp/auxiliary/Random.h" namespace "Aux::Random":
	void _setSeed "Aux::Random::setSeed" (uint64_t, bool)

def setSeed(uint64_t seed, bool useThreadId):
	""" Set the random seed that is used in NetworKit.

	Note that there is a separate random number generator per thread.

	Parameters
	----------
	seed : uint64_t
		The seed
	useThreadId : bool
		If the thread id shall be added to the seed
	"""
	_setSeed(seed, useThreadId)

# Class definitions

## Module: engineering

# TODO: timer

## Module: graph

# DEPRECATED
# TODO: replace with std::pair<double>
cdef extern from "cpp/viz/Point.h" namespace "NetworKit":
	cdef cppclass Point[T]:
		Point()
		Point(T x, T y)
		T& operator[](const index i) except +
		T& at(const index i) except +

cdef extern from "cpp/graph/Graph.h":
	cdef cppclass _Graph "NetworKit::Graph":
		_Graph() except +
		_Graph(count, bool, bool) except +
		_Graph(const _Graph& other) except +
		_Graph(const _Graph& other, bool weighted, bool directed) except +
		void indexEdges() except +
		bool hasEdgeIds() except +
		edgeid edgeId(node, node) except +
		count numberOfNodes() except +
		count numberOfEdges() except +
		index upperNodeIdBound() except +
		index upperEdgeIdBound() except +
		count degree(node u) except +
		count degreeIn(node u) except +
		count degreeOut(node u) except +
		bool isIsolated(node u) except +
		_Graph copyNodes() except +
		node addNode() except +
		void removeNode(node u) except +
		bool hasNode(node u) except +
		void addEdge(node u, node v, edgeweight w) except +
		void setWeight(node u, node v, edgeweight w) except +
		void removeEdge(node u, node v) except +
		void swapEdge(node s1, node t1, node s2, node t2) except +
		void compactEdges() except +
		void sortEdges() except +
		bool hasEdge(node u, node v) except +
		edgeweight weight(node u, node v) except +
		vector[node] nodes() except +
		vector[pair[node, node]] edges() except +
		vector[node] neighbors(node u) except +
		void forEdges[Callback](Callback c) except +
		void forNodes[Callback](Callback c) except +
		void forNodePairs[Callback](Callback c) except +
		void forNodesInRandomOrder[Callback](Callback c) except +
		void forEdgesOf[Callback](node u, Callback c) except +
		void forInEdgesOf[Callback](node u, Callback c) except +
		bool isWeighted() except +
		bool isDirected() except +
		string toString() except +
		string getName() except +
		void setName(string name) except +
		edgeweight totalEdgeWeight() except +
		node randomNode() except +
		node randomNeighbor(node) except +
		pair[node, node] randomEdge() except +
		Point[float] getCoordinate(node v) except +
		void setCoordinate(node v, Point[float] value) except +
		void initCoordinates() except +
		count numberOfSelfLoops() except +
		_Graph toUndirected() except +
		void BFSfromNode "BFSfrom"[Callback] (node r, Callback c) except +
		void BFSfrom[Callback](vector[node] startNodes, Callback c) except +
		void BFSEdgesFrom[Callback](node r, Callback c) except +
		void DFSfrom[Callback](node r, Callback c) except +
		void DFSEdgesFrom[Callback](node r, Callback c) except +

cdef cppclass EdgeCallBackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, node v, edgeweight w, edgeid eid):
		cdef bool error = False
		cdef string message
		try:
			(<object>callback)(u, v, w, eid)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef cppclass NodeCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u):
		cdef bool error = False
		cdef string message
		try:
			(<object>callback)(u)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef cppclass NodeDistCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, count dist):
		cdef bool error = False
		cdef string message
		try:
			(<object>callback)(u, dist)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef cppclass NodePairCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, node v):
		cdef bool error = False
		cdef string message
		try:
			(<object>callback)(u, v)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

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
	cdef _Graph _this

	def __cinit__(self, n=0, bool weighted=False, bool directed=False):
		if isinstance(n, Graph):
			self._this = move(_Graph((<Graph>n)._this, weighted, directed))
		else:
			self._this = move(_Graph(<count>n, weighted, directed))

	cdef setThis(self, _Graph& other):
		swap[_Graph](self._this, other)
		return self

	def __copy__(self):
		"""
		Generates a copy of the graph
		"""
		return Graph().setThis(_Graph(self._this))

	def __deepcopy__(self, memo):
		"""
		Generates a (deep) copy of the graph
		"""
		return Graph().setThis(_Graph(self._this))

	def __str__(self):
		return "NetworKit.Graph(name={0}, n={1}, m={2})".format(self.getName(), self.numberOfNodes(), self.numberOfEdges())


	def copyNodes(self):
		"""
		Copies all nodes to a new graph

		Returns
		-------
		Graph
			Graph with the same nodes (without edges)
		"""
		return Graph().setThis(self._this.copyNodes())

	def indexEdges(self):
		"""
		Assign integer ids to edges.

		"""
		self._this.indexEdges()

	def hasEdgeIds(self):
		"""
		Returns true if edges have been indexed

		Returns
		-------
		bool
			if edges have been indexed
		"""
		return self._this.hasEdgeIds()

	def edgeId(self, node u, node v):
		"""
		Returns
		-------
		edgeid
			id of the edge
		"""
		return self._this.edgeId(u, v)

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

	def upperNodeIdBound(self):
		"""
		Get an upper bound for the node ids in the graph

		Returns
		-------
		count
			An upper bound for the node ids in the graph
		"""
		return self._this.upperNodeIdBound()

	def upperEdgeIdBound(self):
		"""
		Get an upper bound for the edge ids in the graph

		Returns
		-------
		count
			An upper bound for the edge ids in the graph
		"""
		return self._this.upperEdgeIdBound()

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

	def isIsolated(self, u):
		"""
		If the node `u` is isolated

		Parameters
		----------
		u : node
			Node.

		Returns
		-------
		bool
			If the node is isolated
		"""
		return self._this.isIsolated(u)

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

	def hasNode(self, u):
		""" Checks if the Graph has the node `u`, i.e. if `u` hasn't been deleted and is in the range of valid ids.

		Parameters
		----------
		u : node
			Node

		Returns
		-------
		bool
			If the Graph has the node `u`
		"""
		return self._this.hasNode(u)

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

	def setWeight(self, u, v, w):
		""" Set the weight of an edge. If the edge does not exist, it will be inserted.

		Parameters
		----------
		u : node
			Endpoint of edge.
		v : node
			Endpoint of edge.
		w : edgeweight
			Edge weight.
		"""
		self._this.setWeight(u, v, w)

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

	def swapEdge(self, node s1, node t1, node s2, node t2):
		"""
		Changes the edge (s1, t1) into (s1, t2) and the edge (s2, t2) into (s2, t1).

		If there are edge weights or edge ids, they are preserved. Note that no check is performed if the swap is actually possible, i.e. does not generate duplicate edges.

		Parameters
		----------
		s1 : node
			Source node of the first edge
		t1 : node
			Target node of the first edge
		s2 : node
			Source node of the second edge
		t2 : node
			Target node of the second edge
		"""
		self._this.swapEdge(s1, t1, s2, t2)

	def compactEdges(self):
		"""
		Compact the edge storage, this should be called after executing many edge deletions.
		"""
		self._this.compactEdges()

	def sortEdges(self):
		"""
		Sorts the adjacency arrays by node id. While the running time is linear this
		temporarily duplicates the memory.
		"""
		self._this.sortEdges()

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

	def forNodes(self, object callback):
		""" Experimental node iterator interface

		Parameters
		----------
		callback : object
			Any callable object that takes the parameter node
		"""
		cdef NodeCallbackWrapper* wrapper
		try:
			wrapper = new NodeCallbackWrapper(callback)
			self._this.forNodes[NodeCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forNodesInRandomOrder(self, object callback):
		""" Experimental node iterator interface

		Parameters
		----------
		callback : object
			Any callable object that takes the parameter node
		"""
		cdef NodeCallbackWrapper* wrapper
		try:
			wrapper = new NodeCallbackWrapper(callback)
			self._this.forNodesInRandomOrder[NodeCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forNodePairs(self, object callback):
		""" Experimental node pair iterator interface

		Parameters
		----------
		callback : object
			Any callable object that takes the parameters (node, node)
		"""
		cdef NodePairCallbackWrapper* wrapper
		try:
			wrapper = new NodePairCallbackWrapper(callback)
			self._this.forNodePairs[NodePairCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forEdges(self, object callback):
		""" Experimental edge iterator interface

		Parameters
		----------
		callback : object
			Any callable object that takes the parameter (node, node, edgeweight, edgeid)
		"""
		cdef EdgeCallBackWrapper* wrapper
		try:
			wrapper = new EdgeCallBackWrapper(callback)
			self._this.forEdges[EdgeCallBackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forEdgesOf(self, node u, object callback):
		""" Experimental incident (outgoing) edge iterator interface

		Parameters
		----------
		u : node
			The node of which incident edges shall be passed to the callback
		callback : object
			Any callable object that takes the parameter (node, node, edgeweight, edgeid)
		"""
		cdef EdgeCallBackWrapper* wrapper
		try:
			wrapper = new EdgeCallBackWrapper(callback)
			self._this.forEdgesOf[EdgeCallBackWrapper](u, dereference(wrapper))
		finally:
			del wrapper

	def forInEdgesOf(self, node u, object callback):
		""" Experimental incident incoming edge iterator interface

		Parameters
		----------
		u : node
			The node of which incident edges shall be passed to the callback
		callback : object
			Any callable object that takes the parameter (node, node, edgeweight, edgeid)
		"""
		cdef EdgeCallBackWrapper* wrapper
		try:
			wrapper = new EdgeCallBackWrapper(callback)
			self._this.forInEdgesOf[EdgeCallBackWrapper](u, dereference(wrapper))
		finally:
			del wrapper

	def toUndirected(self):
		"""
		Return an undirected version of this graph.

	 	Returns
	 	-------
			undirected graph.
		"""
		return Graph().setThis(self._this.toUndirected())

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

	def getCoordinate(self, v):
		"""
		DEPRECATED: Coordinates should be handled outside the Graph class
		 like general node attributes.

		Get the coordinates of node v.
		Parameters
		----------
		v : node
			Node.

		Returns
		-------
		pair[float, float]
			x and y coordinates of v.
		"""

		return (self._this.getCoordinate(v)[0], self._this.getCoordinate(v)[1])

	def setCoordinate(self, v, value):
		"""
		DEPRECATED: Coordinates should be handled outside the Graph class
		 like general node attributes.

		Set the coordinates of node v.
		Parameters
		----------
		v : node
			Node.
		value : pair[float, float]
			x and y coordinates of v.
		"""
		cdef Point[float] p = Point[float](value[0], value[1])
		self._this.setCoordinate(v, p)

	def initCoordinates(self):
		"""
		DEPRECATED: Coordinates should be handled outside the Graph class
		 like general node attributes.
		"""
		self._this.initCoordinates()

	def numberOfSelfLoops(self):
		""" Get number of self-loops, i.e. edges {v, v}.
		Returns
		-------
		count
			number of self-loops.
		"""
		return self._this.numberOfSelfLoops()

	def BFSfrom(self, start, object callback):
		""" Experimental BFS search interface

		Parameters
		----------
		start: node or list[node]
			One or more start nodes from which the BFS shall be started
		callback : object
			Any callable object that takes the parameter (node, count) (the second parameter is the depth)
		"""
		cdef NodeDistCallbackWrapper *wrapper
		try:
			wrapper = new NodeDistCallbackWrapper(callback)
			try:
				self._this.BFSfromNode[NodeDistCallbackWrapper](<node?>start, dereference(wrapper))
			except TypeError:
				self._this.BFSfrom[NodeDistCallbackWrapper](<vector[node]?>start, dereference(wrapper))
		finally:
			del wrapper

	def BFSEdgesFrom(self, node start, object callback):
		""" Experimental BFS search interface that passes edges that are part of the BFS tree to the callback

		Parameters
		----------
		start: node
			The start node from which the BFS shall be started
		callback : object
			Any callable object that takes the parameter (node, node)
		"""
		cdef NodePairCallbackWrapper *wrapper
		try:
			wrapper = new NodePairCallbackWrapper(callback)
			self._this.BFSEdgesFrom[NodePairCallbackWrapper](start, dereference(wrapper))
		finally:
			del wrapper

	def DFSfrom(self, node start, object callback):
		""" Experimental DFS search interface

		Parameters
		----------
		start: node
			The start node from which the DFS shall be started
		callback : object
			Any callable object that takes the parameter node
		"""
		cdef NodeCallbackWrapper *wrapper
		try:
			wrapper = new NodeCallbackWrapper(callback)
			self._this.DFSfrom[NodeCallbackWrapper](start, dereference(wrapper))
		finally:
			del wrapper

	def DFSEdgesFrom(self, node start, object callback):
		""" Experimental DFS search interface that passes edges that are part of the DFS tree to the callback

		Parameters
		----------
		start: node
			The start node from which the DFS shall be started
		callback : object
			Any callable object that takes the parameter (node, node)
		"""
		cdef NodePairCallbackWrapper *wrapper
		try:
			wrapper = new NodePairCallbackWrapper(callback)
			self._this.DFSEdgesFrom(start, dereference(wrapper))
		finally:
			del wrapper

# TODO: expose all methods

cdef extern from "cpp/graph/BFS.h":
	cdef cppclass _BFS "NetworKit::BFS":
		_BFS(_Graph G, node source, bool storePaths, bool storeStack) except +
		void run() nogil except +
		void run(node t) nogil except +
		vector[edgeweight] getDistances() except +
		vector[node] getPath(node t) except +

cdef class BFS:
	""" Simple breadth-first search on a Graph from a given source

	BFS(G, source, [storePaths], [storeStack])

	Create BFS for `G` and source node `source`.

	Parameters
	----------
	G : Graph
		The graph.
	source : node
		The source node of the breadth-first search.
	storePaths : bool
		store paths and number of paths?

	"""
	cdef _BFS* _this
	cdef Graph _G

	def __cinit__(self, Graph G, source, storePaths=True, storeStack=False):
		self._G = G
		self._this = new _BFS(G._this, source, storePaths, storeStack)

	def __dealloc__(self):
		del self._this

	def run(self, t = None):
		"""
		Breadth-first search from source.

		Returns
		-------
		vector
			Vector of unweighted distances from source node, i.e. the
	 		length (number of edges) of the shortest path from source to any other node.
		"""
		cdef node ct
		if t == None:
			with nogil:
				self._this.run()
		else:
			ct = <node>t
			with nogil:
				self._this.run(ct)
		return self

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


cdef extern from "cpp/graph/DynBFS.h":
	cdef cppclass _DynBFS "NetworKit::DynBFS":
		_DynBFS(_Graph G, node source) except +
		void run() nogil except +
		vector[edgeweight] getDistances() except +
		vector[node] getPath(node t) except +
		void update(vector[_GraphEvent]) except +

cdef class DynBFS:
	""" Dynamic version of BFS.

	DynBFS(G, source)

	Create DynBFS for `G` and source node `source`.

	Parameters
	----------
	G : Graph
		The graph.
	source : node
		The source node of the breadth-first search.
	storeStack : bool
		maintain a stack of nodes in order of decreasing distance?
	"""
	cdef _DynBFS* _this
	cdef Graph _G

	def __cinit__(self, Graph G, source):
		self._G = G
		self._this = new _DynBFS(G._this, source)

	# this is necessary so that the C++ object gets properly garbage collected
	def __dealloc__(self):
		del self._this

	def run(self):
		"""
		Breadth-first search from source.

		Returns
		-------
		vector
			Vector of unweighted distances from source node, i.e. the
			length (number of edges) of the shortest path from source to any other node.
		"""
		with nogil:
			self._this.run()
		return self

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

	def update(self, batch):
		""" Updates shortest paths with the batch `batch` of edge insertions.

		Parameters
		----------
		batch : list of GraphEvent.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.update(_batch)



cdef extern from "cpp/graph/Dijkstra.h":
	cdef cppclass _Dijkstra "NetworKit::Dijkstra":
		_Dijkstra(_Graph G, node source, bool storePaths, bool storeStack) except +
		void run() nogil except +
		void run(node t) nogil except +
		vector[edgeweight] getDistances() except +
		vector[node] getPath(node t) except +

cdef class Dijkstra:
	""" Dijkstra's SSSP algorithm.
	Returns list of weighted distances from node source, i.e. the length of the shortest path from source to
	any other node.

    Dijkstra(G, source, [storePaths], [storeStack])

    Creates Dijkstra for `G` and source node `source`.

    Parameters
	----------
	G : Graph
		The graph.
	source : node
		The source node.
	storePaths : bool
		store paths and number of paths?
	storeStack : bool
		maintain a stack of nodes in order of decreasing distance?
    """
	cdef _Dijkstra* _this
	cdef Graph _G

	def __cinit__(self, Graph G, source, storePaths=True, storeStack=False):
		self._G = G
		self._this = new _Dijkstra(G._this, source, storePaths, storeStack)

	def __dealloc__(self):
		del self._this

	def run(self, t = None):
		"""
		Breadth-first search from source.

		Returns
		-------
		vector
			Vector of unweighted distances from source node, i.e. the
	 		length (number of edges) of the shortest path from source to any other node.
		"""
		cdef node ct
		if t == None:
			with nogil:
				self._this.run()
		else:
			ct = <node>t
			with nogil:
				self._this.run(ct)
		return self

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


cdef extern from "cpp/graph/DynDijkstra.h":
	cdef cppclass _DynDijkstra "NetworKit::DynDijkstra":
		_DynDijkstra(_Graph G, node source) except +
		void run() nogil except +
		vector[edgeweight] getDistances() except +
		vector[node] getPath(node t) except +
		void update(vector[_GraphEvent]) except +

cdef class DynDijkstra:
	""" Dynamic version of Dijkstra.

	DynDijkstra(G, source)

	Create DynDijkstra for `G` and source node `source`.

	Parameters
	----------
	G : Graph
		The graph.
	source : node
		The source node of the breadth-first search.

	"""
	cdef _DynDijkstra* _this
	cdef Graph _G

	def __cinit__(self, Graph G, source):
		self._G = G
		self._this = new _DynDijkstra(G._this, source)

	# this is necessary so that the C++ object gets properly garbage collected
	def __dealloc__(self):
		del self._this

	def init(self):
		"""
		SSSP search from source.

		Returns
		-------
		vector
			Vector of distances from source node, i.e. the length of the
			shortest path from source to any other node.
		"""
		with nogil:
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

	def update(self, batch):
		""" Updates shortest paths with the batch `batch` of edge insertions.

		Parameters
		----------
		batch : list of GraphEvent.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.update(_batch)


cdef extern from "cpp/graph/Subgraph.h" namespace "NetworKit::Subgraph":
		_Graph _SubGraphFromNodes "NetworKit::Subgraph::fromNodes"(_Graph G, unordered_set[node] nodes)  except +

cdef class Subgraph:
	""" Methods for creating subgraphs """

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
		return Graph().setThis(_SubGraphFromNodes(G._this, nnodes))


cdef extern from "cpp/graph/SpanningForest.h":
	cdef cppclass _SpanningForest "NetworKit::SpanningForest":
		_SpanningForest(_Graph) except +
		_Graph generate() except +

cdef class SpanningForest:
	""" Generates a spanning forest for a given graph

		Parameters
		----------
		G : Graph
			The graph.
		nodes : list
			A subset of nodes of `G` which induce the subgraph.
	"""
	cdef _SpanningForest* _this
	cdef Graph _G

	def __cinit__(self, Graph G not None):
		self._G = G
		self._this = new _SpanningForest(G._this)

	def __dealloc__(self):
		del self._this

	def generate(self):
		return Graph().setThis(self._this.generate());




cdef extern from "cpp/independentset/Luby.h":
	cdef cppclass _Luby "NetworKit::Luby":
		_Luby() except +
		vector[bool] run(_Graph G) except +
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
		return self._this.run(G._this)
		# TODO: return self

	def toString(self):
		""" Get string representation of the algorithm.

		Returns
		-------
		string
			The string representation of the algorithm.
		"""
		return self._this.toString().decode("utf-8")


# Module: generators

# cdef extern from "cpp/graph/GraphGenerator.h":
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

cdef extern from "cpp/generators/BarabasiAlbertGenerator.h":
	cdef cppclass _BarabasiAlbertGenerator "NetworKit::BarabasiAlbertGenerator":
		_BarabasiAlbertGenerator() except +
		_BarabasiAlbertGenerator(count k, count nMax, count n0) except +
		#_Graph* _generate()
		_Graph generate() except +

cdef class BarabasiAlbertGenerator:
	""" Generates a scale-free graph using the Barabasi-Albert preferential attachment model.

	Parameters
	----------
	k : count
		number of edges that come with a new node
	nMax : count
		maximum number of nodes produced
	n0 : count
		number of starting nodes
	 """
	cdef _BarabasiAlbertGenerator _this

	def __cinit__(self, k, nMax, n0):
		self._this = _BarabasiAlbertGenerator(k, nMax, n0)

	def generate(self):
		return Graph().setThis(self._this.generate());


cdef extern from "cpp/generators/PubWebGenerator.h":
	cdef cppclass _PubWebGenerator "NetworKit::PubWebGenerator":
		_PubWebGenerator(count numNodes, count numberOfDenseAreas, float neighborhoodRadius, count maxNumberOfNeighbors) except +
		_Graph generate() except +

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
		self._this = new _PubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	def __dealloc__(self):
		del self._this

	def generate(self):
		return Graph(0).setThis(self._this.generate())


cdef extern from "cpp/generators/ErdosRenyiGenerator.h":
	cdef cppclass _ErdosRenyiGenerator "NetworKit::ErdosRenyiGenerator":
		_ErdosRenyiGenerator(count nNodes, double prob, bool directed) except +
		_Graph generate() except +

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
	directed : bool
		Generates a directed
	"""

	cdef _ErdosRenyiGenerator* _this

	def __cinit__(self, nNodes, prob, directed=False):
		self._this = new _ErdosRenyiGenerator(nNodes, prob, directed)

	def __dealloc__(self):
		del self._this

	def generate(self):
		return Graph(0).setThis(self._this.generate())

cdef extern from "cpp/generators/DorogovtsevMendesGenerator.h":
	cdef cppclass _DorogovtsevMendesGenerator "NetworKit::DorogovtsevMendesGenerator":
		_DorogovtsevMendesGenerator(count nNodes) except +
		_Graph generate() except +

cdef class DorogovtsevMendesGenerator:
	""" Generates a graph according to the Dorogovtsev-Mendes model.

 	DorogovtsevMendesGenerator(nNodes)

 	Constructs the generator class.

	Parameters
	----------
	nNodes : count
		Number of nodes in the target graph.
	"""

	cdef _DorogovtsevMendesGenerator* _this

	def __cinit__(self, nNodes):
		self._this = new _DorogovtsevMendesGenerator(nNodes)

	def __dealloc__(self):
		del self._this

	def generate(self):
		""" Generates a random graph according to the Dorogovtsev-Mendes model.

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph(0).setThis(self._this.generate())


cdef extern from "cpp/generators/RegularRingLatticeGenerator.h":
	cdef cppclass _RegularRingLatticeGenerator "NetworKit::RegularRingLatticeGenerator":
		_RegularRingLatticeGenerator(count nNodes, count nNeighbors) except +
		_Graph generate() except +

cdef class RegularRingLatticeGenerator:
	"""
	Constructs a regular ring lattice.

	RegularRingLatticeGenerator(count nNodes, count nNeighbors)

	Constructs the generator.

	Parameters
	----------
	nNodes : number of nodes in the target graph.
	nNeighbors : number of neighbors on each side of a node
	"""

	cdef _RegularRingLatticeGenerator* _this

	def __cinit__(self, nNodes, nNeighbors):
		self._this = new _RegularRingLatticeGenerator(nNodes, nNeighbors)

	def __dealloc__(self):
		del self._this

	def generate(self):
		""" Generates a rgular ring lattice.

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph(0).setThis(self._this.generate())


cdef extern from "cpp/generators/WattsStrogatzGenerator.h":
	cdef cppclass _WattsStrogatzGenerator "NetworKit::WattsStrogatzGenerator":
		_WattsStrogatzGenerator(count nNodes, count nNeighbors, double p) except +
		_Graph generate() except +

cdef class WattsStrogatzGenerator:
	""" Generates a graph according to the Watts-Strogatz model.

	First, a regular ring lattice is generated. Then edges are rewired
		with a given probability.

	WattsStrogatzGenerator(count nNodes, count nNeighbors, double p)

	Constructs the generator.

	Parameters
	----------
	nNodes : Number of nodes in the target graph.
	nNeighbors : number of neighbors on each side of a node
	p : rewiring probability
	"""

	cdef _WattsStrogatzGenerator* _this

	def __dealloc__(self):
		del self._this

	def __cinit__(self, nNodes, nNeighbors, p):
		self._this = new _WattsStrogatzGenerator(nNodes, nNeighbors, p)

	def generate(self):
		""" Generates a random graph according to the Watts-Strogatz model.

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph(0).setThis(self._this.generate())


cdef extern from "cpp/generators/ClusteredRandomGraphGenerator.h":
	cdef cppclass _ClusteredRandomGraphGenerator "NetworKit::ClusteredRandomGraphGenerator":
		_ClusteredRandomGraphGenerator(count, count, double, double) except +
		_Graph generate() except +
		_Partition getCommunities() except +

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

	def __dealloc__(self):
		del self._this

	def generate(self):
		""" Generates a clustered random graph with the properties given in the constructor.

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph(0).setThis(self._this.generate())

	def getCommunities(self):
		""" Returns the generated ground truth clustering.

		Returns
		-------
		Partition
			The generated ground truth clustering.
		"""
		return Partition().setThis(self._this.getCommunities())


cdef extern from "cpp/generators/ChungLuGenerator.h":
	cdef cppclass _ChungLuGenerator "NetworKit::ChungLuGenerator":
		_ChungLuGenerator(vector[count] degreeSequence) except +
		_Graph generate() except +

cdef class ChungLuGenerator:
	"""
		Given an arbitrary degree sequence, the Chung-Lu generative model
		will produce a random graph with the same expected degree sequence.

		see Chung, Lu: The average distances in random graphs with given expected degrees
		and Chung, Lu: Connected Components in Random Graphs with Given Expected Degree Sequences.
		Aiello, Chung, Lu: A Random Graph Model for Massive Graphs describes a different generative model
		which is basically asymptotically equivalent but produces multi-graphs.
	"""

	cdef _ChungLuGenerator* _this

	def __cinit__(self, vector[count] degreeSequence):
		self._this = new _ChungLuGenerator(degreeSequence)

	def __dealloc__(self):
		del self._this

	def generate(self):
		""" Generates graph with expected degree sequence seq.

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph(0).setThis(self._this.generate())


cdef extern from "cpp/generators/HavelHakimiGenerator.h":
	cdef cppclass _HavelHakimiGenerator "NetworKit::HavelHakimiGenerator":
		_HavelHakimiGenerator(vector[count] degreeSequence, bool ignoreIfRealizable) except +
		_Graph generate() except +
		bool isRealizable() except +
		bool getRealizable() except +

cdef class HavelHakimiGenerator:
	""" Havel-Hakimi algorithm for generating a graph according to a given degree sequence.

		The sequence, if it is realizable, is reconstructed exactly. The resulting graph usually
		has a high clustering coefficient. Construction runs in linear time O(m).

		If the sequence is not realizable, depending on the parameter ignoreIfRealizable, either
		an exception is thrown during generation or the graph is generated with a modified degree
		sequence, i.e. not all nodes might have as many neighbors as requested.

		HavelHakimiGenerator(sequence, ignoreIfRealizable=True)

		Parameters
		----------
		sequence : vector
			Degree sequence to realize. Must be non-increasing.
		ignoreIfRealizable : bool, optional
			If true, generate the graph even if the degree sequence is not realizable. Some nodes may get lower degrees than requested in the sequence.
	"""

	cdef _HavelHakimiGenerator* _this


	def __cinit__(self, vector[count] degreeSequence, ignoreIfRealizable=True):
		self._this = new _HavelHakimiGenerator(degreeSequence, ignoreIfRealizable)

	def __dealloc__(self):
		del self._this

	def isRealizable(self):
		return self._this.isRealizable()

	def getRealizable(self):
		return self._this.getRealizable();

	def generate(self):
		""" Generates degree sequence seq (if it is realizable).

		Returns
		-------
		Graph
			Graph with degree sequence seq or modified sequence if ignoreIfRealizable is true and the sequence is not realizable.
		"""
		return Graph(0).setThis(self._this.generate())

cdef extern from "cpp/generators/ConfigurationModelGenerator.h":
	cdef cppclass _ConfigurationModelGenerator "NetworKit::ConfigurationModelGenerator":
		_ConfigurationModelGenerator(vector[count] degreeSequence, bool ignoreIfRealizable) except +
		_Graph generate() except +
		bool isRealizable() except +
		bool getRealizable() except +

cdef class ConfigurationModelGenerator:
	"""
	Configuration model graph generator for generating a random simple graph with exactly the given degree sequence.

	This implementation is based on the paper
	"Random generation of large connected simple graphs with prescribed degree distribution" by Fabien Viger and Matthieu Latapy,
	available at http://www-rp.lip6.fr/~latapy/FV/generation.html, however without preserving connectivity (this could later be added as
	optional feature).

	The Havel-Hakami generator is used for the initial graph generation, then the Markov-Chain Monte-Carlo algorithm as described and
	implemented by Fabien Viger and Matthieu Latapy but without the steps for ensuring connectivity is executed. This should lead to a
	graph that is drawn uniformly at random from all graphs with the given degree sequence.

	Note that at most 10 times the number of edges edge swaps are performed (same number as in the abovementioned implementation) and
	in order to limit the running time, at most 200 times as many attempts to perform an edge swap are made (as certain degree distributions
	do not allow edge swaps at all).

	Parameters
	----------
	degreeSequence : vector[count]
		The degree sequence that shall be generated
	ignoreIfRealizable : bool, optional
		If true, generate the graph even if the degree sequence is not realizable. Some nodes may get lower degrees than requested in the sequence.
	"""
	cdef _ConfigurationModelGenerator *_this

	def __cinit__(self, vector[count] degreeSequence, bool ignoreIfRealizable = False):
		self._this = new _ConfigurationModelGenerator(degreeSequence, ignoreIfRealizable)

	def __dealloc__(self):
		del self._this

	def isRealizable(self):
		return self._this.isRealizable()

	def getRealizable(self):
		return self._this.getRealizable()

	def generate(self):
		"""
		Generate a graph according to the configuration model.

		Issues a INFO log message if the wanted number of edge swaps cannot be performed because of the limit of attempts (see in the description of the class for details).

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph().setThis(self._this.generate())


cdef extern from "cpp/generators/RmatGenerator.h":
	cdef cppclass _RmatGenerator "NetworKit::RmatGenerator":
		_RmatGenerator(count scale, count edgeFactor, double a, double b, double c, double d) except +
		_Graph generate() except +

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

	def __dealloc__(self):
		del self._this

	def generate(self):
		""" Graph to be generated according to parameters specified in constructor.

		Returns
		-------
		Graph
			The generated graph.
		"""
		return Graph(0).setThis(self._this.generate())

cdef extern from "cpp/generators/PowerlawDegreeSequence.h":
	cdef cppclass _PowerlawDegreeSequence "NetworKit::PowerlawDegreeSequence":
		_PowerlawDegreeSequence(count minDeg, count maxDeg, double gamma) except +
		void setMinimumFromAverageDegree(double avgDeg) nogil except +
		double getExpectedAverageDegree() except +
		count getMinimumDegree() const
		void run() nogil except +
		vector[count] getDegreeSequence(count numNodes) except +
		count getDegree() except +

cdef class PowerlawDegreeSequence:
	"""
	Generates a powerlaw degree sequence with the given minimum and maximum degree, the powerlaw exponent gamma

	Parameters
	----------
	minDeg : count
		The minium degree
	maxDeg : count
		The maximum degree
	gamma : double
		The powerlaw exponent
	"""
	cdef _PowerlawDegreeSequence *_this

	def __cinit__(self, count minDeg, count maxDeg, double gamma):
		self._this = new _PowerlawDegreeSequence(minDeg, maxDeg, gamma)

	def __dealloc__(self):
		del self._this

	def setMinimumFromAverageDegree(self, double avgDeg):
		"""
		Tries to set the minimum degree such that the specified average degree is expected.

		Parameters
		----------
		avgDeg : double
			The average degree that shall be approximated
		"""
		with nogil:
			self._this.setMinimumFromAverageDegree(avgDeg)
		return self

	def getExpectedAverageDegree(self):
		"""
		Returns the expected average degree. Note: run needs to be called first.

		Returns
		-------
		double
			The expected average degree.
		"""
		return self._this.getExpectedAverageDegree()

	def getMinimumDegree(self):
		"""
		Returns the minimum degree.

		Returns
		-------
		count
			The minimum degree
		"""
		return self._this.getMinimumDegree()

	def run(self):
		"""
		Executes the generation of the probability distribution.
		"""
		with nogil:
			self._this.run()
		return self

	def getDegreeSequence(self, count numNodes):
		"""
		Returns a degree sequence with even degree sum.

		Parameters
		----------
		numNodes : count
			The number of nodes/degrees that shall be returned

		Returns
		-------
		vector[count]
			The generated degree sequence
		"""
		return self._this.getDegreeSequence(numNodes)

	def getDegree(self):
		"""
		Returns a degree drawn at random with a power law distribution

		Returns
		-------
		count
			The generated random degree
		"""
		return self._this.getDegree()

# Module: graphio

cdef extern from "cpp/io/GraphReader.h":
	cdef cppclass _GraphReader "NetworKit::GraphReader":
		_GraphReader() nogil except +
		_Graph read(string path) nogil except +

cdef class GraphReader:
	""" Abstract base class for graph readers"""

	cdef _GraphReader* _this

	def __init__(self, *args, **kwargs):
		if type(self) == GraphReader:
			raise RuntimeError("Error, you may not use GraphReader directly, use a sub-class instead")

	def __cinit__(self, *args, **kwargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

	def read(self, path):
		cdef string cpath = stdstring(path)
		cdef _Graph result

		with nogil:
			result = move(self._this.read(cpath)) # extra move in order to avoid copying the internal variable that is used by Cython
		return Graph(0).setThis(result)

cdef extern from "cpp/io/METISGraphReader.h":
	cdef cppclass _METISGraphReader "NetworKit::METISGraphReader" (_GraphReader):
		_METISGraphReader() nogil except +

cdef class METISGraphReader(GraphReader):
	""" Reads the METIS adjacency file format [1]. If the Fast reader fails,
		use readGraph(path, graphio.formats.metis) as an alternative.
		[1]: http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
	"""
	def __cinit__(self):
		self._this = new _METISGraphReader()

cdef extern from "cpp/io/GraphToolBinaryReader.h":
	cdef cppclass _GraphToolBinaryReader "NetworKit::GraphToolBinaryReader" (_GraphReader):
		_GraphToolBinaryReader() except +

cdef class GraphToolBinaryReader(GraphReader):
	""" Reads the binary file format defined by graph-tool[1].
		[1]: http://graph-tool.skewed.de/static/doc/gt_format.html
	"""
	def __cinit__(self):
		self._this = new _GraphToolBinaryReader()


cdef extern from "cpp/io/EdgeListReader.h":
	cdef cppclass _EdgeListReader "NetworKit::EdgeListReader"(_GraphReader):
		_EdgeListReader() except +
		_EdgeListReader(char separator, node firstNode, string commentPrefix, bool continuous, bool directed)
		unordered_map[node,node] getNodeMap() except +


cdef class EdgeListReader(GraphReader):
	""" Reads a file in an edge list format.
		TODO: docstring
	"""
	def __cinit__(self, separator, firstNode, commentPrefix="#", continuous=True, directed=False):
		self._this = new _EdgeListReader(stdstring(separator)[0], firstNode, stdstring(commentPrefix), continuous, directed)

	def getNodeMap(self):
		cdef unordered_map[node,node] cResult = (<_EdgeListReader*>(self._this)).getNodeMap()
		result = []
		for elem in cResult:
			result.append((elem.first,elem.second))
		return result

cdef extern from "cpp/io/KONECTGraphReader.h":
	cdef cppclass _KONECTGraphReader "NetworKit::KONECTGraphReader"(_GraphReader):
		_KONECTGraphReader() except +
		_KONECTGraphReader(char separator, bool ignoreLoops)

cdef class KONECTGraphReader(GraphReader):
	""" Reader for the KONECT graph format, which is described in detail on the KONECT website[1].

		[1]: http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
	"""
	def __cinit__(self, separator, ignoreLoops = False):
		self._this = new _KONECTGraphReader(stdstring(separator)[0], ignoreLoops)

cdef extern from "cpp/io/GMLGraphReader.h":
	cdef cppclass _GMLGraphReader "NetworKit::GMLGraphReader"(_GraphReader):
		_GMLGraphReader() except +

cdef class GMLGraphReader(GraphReader):
	""" Reader for the GML graph format, which is documented here [1].

		[1]: http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf
 	"""
	def __cinit__(self):
		self._this = new _GMLGraphReader()

cdef extern from "cpp/io/METISGraphWriter.h":
	cdef cppclass _METISGraphWriter "NetworKit::METISGraphWriter":
		_METISGraphWriter() except +
		void write(_Graph G, string path) nogil except +


cdef class METISGraphWriter:
	""" Writes graphs in the METIS format"""
	cdef _METISGraphWriter _this

	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(G._this, cpath)

cdef extern from "cpp/io/GraphToolBinaryWriter.h":
	cdef cppclass _GraphToolBinaryWriter "NetworKit::GraphToolBinaryWriter":
		_GraphToolBinaryWriter() except +
		void write(_Graph G, string path) nogil except +


cdef class GraphToolBinaryWriter:
	""" Reads the binary file format defined by graph-tool[1].
		[1]: http://graph-tool.skewed.de/static/doc/gt_format.html
	"""
	cdef _GraphToolBinaryWriter _this

	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(G._this, cpath)


cdef extern from "cpp/io/DotGraphWriter.h":
	cdef cppclass _DotGraphWriter "NetworKit::DotGraphWriter":
		_DotGraphWriter() except +
		void write(_Graph G, string path) nogil except +


cdef class DotGraphWriter:
	""" Writes graphs in the .dot/GraphViz format"""
	cdef _DotGraphWriter _this

	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(G._this, cpath)


#cdef extern from "cpp/io/VNAGraphWriter.h":
#	cdef cppclass _VNAGraphWriter "NetworKit::VNAGraphWriter":
#		_VNAGraphWriter() except +
#		void write(_Graph G, string path) except +


#cdef class VNAGraphWriter:
#	""" Writes graphs in the VNA format. The VNA format is commonly used by Netdraw, and is very similar to Pajek format.
#	It defines nodes and edges (ties), and supports attributes. Each section of the file is separated by an asterisk. """
#	cdef _VNAGraphWriter _this

#	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
#		self._this.write(G._this, stdstring(path))


cdef extern from "cpp/io/GMLGraphWriter.h":
	cdef cppclass _GMLGraphWriter "NetworKit::GMLGraphWriter":
		_GMLGraphWriter() except +
		void write(_Graph G, string path) nogil except +


cdef class GMLGraphWriter:
	""" Writes a graph and its coordinates as a GML file.[1]
		[1] http://svn.bigcat.unimaas.nl/pvplugins/GML/trunk/docs/gml-technical-report.pdf """
	cdef _GMLGraphWriter _this

	def write(self, Graph G not None, path):
		 # string needs to be converted to bytes, which are coerced to std::string
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(G._this, cpath)


cdef extern from "cpp/io/EdgeListWriter.h":
	cdef cppclass _EdgeListWriter "NetworKit::EdgeListWriter":
		_EdgeListWriter() except +
		_EdgeListWriter(char separator, node firstNode) except +
		void write(_Graph G, string path) nogil except +

cdef class EdgeListWriter:
	""" Reads and writes graphs in various edge list formats. The constructor takes a
		seperator char and the ID of the first node as paraneters."""

	cdef _EdgeListWriter _this

	def __cinit__(self, separator, firstNode):
		cdef char sep = stdstring(separator)[0]
		self._this = _EdgeListWriter(sep, firstNode)

	def write(self, Graph G not None, path):
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(G._this, cpath)



cdef extern from "cpp/io/LineFileReader.h":
	cdef cppclass _LineFileReader "NetworKit::LineFileReader":
		_LineFileReader() except +
		vector[string] read(string path)


cdef class LineFileReader:
	""" Reads a file and puts each line in a list of strings """
	cdef _LineFileReader _this

	def read(self, path):
		return self._this.read(stdstring(path))


cdef extern from "cpp/io/SNAPGraphWriter.h":
	cdef cppclass _SNAPGraphWriter "NetworKit::SNAPGraphWriter":
		_SNAPGraphWriter() except +
		void write(_Graph G, string path) nogil except +

cdef class SNAPGraphWriter:
	""" Writes graphs in a format suitable for the Georgia Tech SNAP software [1]
		[1]: http://snap-graph.sourceforge.net/
	"""
	cdef _SNAPGraphWriter _this

	def write(self, Graph G, path):
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(G._this, cpath)


cdef extern from "cpp/io/SNAPGraphReader.h":
	cdef cppclass _SNAPGraphReader "NetworKit::SNAPGraphReader"(_GraphReader):
		_SNAPGraphReader() except +
		unordered_map[node,node] getNodeIdMap() except +

cdef class SNAPGraphReader(GraphReader):
	""" Reads a graph from the SNAP graph data collection [1] (currently experimental)
		[1]: http://snap.stanford.edu/data/index.html
	"""
	def __cinit__(self):
		self._this = new _SNAPGraphReader()

	def getNodeIdMap(self):
		cdef unordered_map[node,node] cResult = (<_SNAPGraphReader*>(self._this)).getNodeIdMap()
		result = []
		for elem in cResult:
			result.append((elem.first,elem.second))
		return result


cdef extern from "cpp/io/PartitionReader.h":
	cdef cppclass _PartitionReader "NetworKit::PartitionReader":
		_PartitionReader() except +
		_Partition read(string path) except +


cdef class PartitionReader:
	""" Reads a partition from a file.
		File format: line i contains subset id of element i.
	 """
	cdef _PartitionReader _this

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))


cdef extern from "cpp/io/PartitionWriter.h":
	cdef cppclass _PartitionWriter "NetworKit::PartitionWriter":
		_PartitionWriter() except +
		void write(_Partition, string path) nogil except +


cdef class PartitionWriter:
	""" Writes a partition to a file.
		File format: line i contains subset id of element i.
	 """
	cdef _PartitionWriter _this

	def write(self, Partition zeta, path):
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(zeta._this, cpath)


cdef extern from "cpp/io/EdgeListPartitionReader.h":
	cdef cppclass _EdgeListPartitionReader "NetworKit::EdgeListPartitionReader":
		_EdgeListPartitionReader() except +
		_EdgeListPartitionReader(node firstNode, char sepChar) except +
		_Partition read(string path) except +


cdef class EdgeListPartitionReader:
	""" Reads a partition from an edge list type of file
	 """
	cdef _EdgeListPartitionReader _this

	def __cinit__(self, node firstNode=1, sepChar = '\t'):
		self._this = _EdgeListPartitionReader(firstNode, stdstring(sepChar)[0])

	def read(self, path):
		return Partition().setThis(self._this.read(stdstring(path)))

cdef extern from "cpp/io/SNAPEdgeListPartitionReader.h":
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
		return Cover().setThis(self._this.read(stdstring(path), cNodeMap, G._this))

#	def readWithInfo(self,path,nNodes):
#		return Partition().setThis(self._this.readWithInfo(stdstring(path),nNodes))

#not existing yet, maybe in the future?
#cdef extern from "cpp/io/EdgeListPartitionWriter.h":
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

cdef extern from "cpp/io/CoverReader.h":
	cdef cppclass _CoverReader "NetworKit::CoverReader":
		_CoverReader() except +
		_Cover read(string path,_Graph G) except +

cdef class CoverReader:
	""" Reads a cover from a file
		File format: each line contains the space-separated node ids of a community
	 """
	cdef _CoverReader _this

	def read(self, path, Graph G):
		return Cover().setThis(self._this.read(stdstring(path), G._this))

cdef extern from "cpp/io/CoverWriter.h":
	cdef cppclass _CoverWriter "NetworKit::CoverWriter":
		_CoverWriter() except +
		void write(_Cover, string path) nogil except +


cdef class CoverWriter:
	""" Writes a partition to a file.
		File format: each line contains the space-separated node ids of a community
	 """
	cdef _CoverWriter _this

	def write(self, Cover zeta, path):
		cdef string cpath = stdstring(path)
		with nogil:
			self._this.write(zeta._this, cpath)

cdef extern from "cpp/io/EdgeListCoverReader.h":
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
		return Cover().setThis(self._this.read(stdstring(path), G._this))

# Parameters

cdef extern from "cpp/base/Parameters.h":
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
cdef extern from "cpp/structures/Partition.h":
	cdef cppclass _Partition "NetworKit::Partition":
		_Partition() except +
		_Partition(index) except +
		_Partition(_Partition) except +
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
		index operator[](index) except +


cdef class Partition:
	""" Implements a partition of a set, i.e. a subdivision of the
 		set into disjoint subsets.

 		Partition(z=0)

 		Create a new partition data structure for `z` elements.

		Parameters
		----------
		size : index, optional
			Maximum index of an element. Default is 0.
	"""
	cdef _Partition _this

	def __cinit__(self, index size=0):
		self._this = move(_Partition(size))

	def __len__(self):
		"""
		Returns
		-------
		count
			Number of elements in the partition.
		"""
		return self._this.numberOfElements()

	def __getitem__(self, index e):
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

	def __setitem__(self, index e, index s):
		""" Set the set (id) in which the element `e` is contained.

		Parameters
		----------
		e : index
			Index of the element
		s : index
			Index of the subset
		"""
		self._this.addToSubset(s, e)

	def __copy__(self):
		"""
		Generates a copy of the partition
		"""
		return Partition().setThis(_Partition(self._this))

	def __deepcopy__(self):
		"""
		Generates a copy of the partition
		"""
		return Partition().setThis(_Partition(self._this))

	cdef setThis(self,  _Partition& other):
		swap[_Partition](self._this,  other)
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
		return self._this.extend()

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

	def __getitem__(self, index):
		return self._this[index]


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


cdef extern from "cpp/structures/Cover.h":
	cdef cppclass _Cover "NetworKit::Cover":
		_Cover() except +
		_Cover(_Partition p) except +
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

	def __cinit__(self, Partition p = None):
		if p is not None:
			self._this = move(_Cover(p._this))

	cdef setThis(self, _Cover& other):
		swap[_Cover](self._this, other)
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

# Fused type for methods that accept both a partition and a cover
ctypedef fused PartitionCover:
	Partition
	Cover

cdef extern from "cpp/community/ClusteringGenerator.h":
	cdef cppclass _ClusteringGenerator "NetworKit::ClusteringGenerator":
		_ClusteringGenerator() except +
		_Partition makeSingletonClustering(_Graph G) except +
		_Partition makeOneClustering(_Graph G) except +
		_Partition makeRandomClustering(_Graph G, count k) except +
		_Partition makeContinuousBalancedClustering(_Graph G, count k) except +
		_Partition makeNoncontinuousBalancedClustering(_Graph G, count k) except +

cdef class ClusteringGenerator:
	""" Generators for various clusterings """
	cdef _ClusteringGenerator _this
	def makeSingletonClustering(self, Graph G):
		"""  Generate a clustering where each node has its own cluster

		Parameters
		----------
		G: Graph
			The graph for which the clustering shall be generated

		Returns
		-------
		Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeSingletonClustering(G._this))
	def makeOneClustering(self, Graph G):
		"""  Generate a clustering with one cluster consisting of all nodes

		Parameters
		----------
		G: Graph
			The graph for which the clustering shall be generated

		Returns
		-------
		Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeOneClustering(G._this))
	def makeRandomClustering(self, Graph G, count k):
		"""  Generate a clustering with `k` clusters to which nodes are assigned randomly

		Parameters
		----------
		G: Graph
			The graph for which the clustering shall be generated
		k: count
			The number of clusters that shall be generated

		Returns
		-------
		Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeRandomClustering(G._this, k))
	def makeContinuousBalancedClustering(self, Graph G, count k):
		"""  Generate a clustering with `k` clusters to which nodes are assigned in continuous blocks

		Parameters
		----------
		G: Graph
			The graph for which the clustering shall be generated
		k: count
			The number of clusters that shall be generated

		Returns
		-------
		Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeContinuousBalancedClustering(G._this, k))
	def makeNoncontinuousBalancedClustering(self, Graph G, count k):
		"""  Generate a clustering with `k` clusters, the ith node is assigned to cluster i % k. This means that
		for k**2 nodes, this clustering is complementary to the continuous clustering in the sense that no pair
		of nodes that is in the same cluster in one of the clusterings is in the same cluster in the other clustering.

		Parameters
		----------
		G: Graph
			The graph for which the clustering shall be generated
		k: count
			The number of clusters that shall be generated

		Returns
		-------
		Partition
			The generated partition
		"""
		return Partition().setThis(self._this.makeNoncontinuousBalancedClustering(G._this, k))

cdef extern from "cpp/community/GraphClusteringTools.h" namespace "NetworKit::GraphClusteringTools":
	float getImbalance(_Partition zeta) except +
	_Graph communicationGraph(_Graph graph, _Partition zeta) except +
	count weightedDegreeWithCluster(_Graph graph, _Partition zeta, node u, index cid)
	bool isProperClustering(_Graph G, _Partition zeta)
	bool isSingletonClustering(_Graph G, _Partition zeta)
	bool isOneClustering(_Graph G, _Partition zeta)
	bool equalClusterings(_Partition zeta, _Partition eta, _Graph G)

cdef class GraphClusteringTools:
	@staticmethod
	def getImbalance(Partition zeta):
		return getImbalance(zeta._this)
	@staticmethod
	def communicationGraph(Graph graph, Partition zeta):
		return Graph().setThis(communicationGraph(graph._this, zeta._this))
	@staticmethod
	def weightedDegreeWithCluster(Graph graph, Partition zeta, node u, index cid):
		return weightedDegreeWithCluster(graph._this, zeta._this, u, cid)
	@staticmethod
	def isProperClustering(Graph G, Partition zeta):
		return isProperClustering(G._this, zeta._this)
	@staticmethod
	def isSingletonClustering(Graph G, Partition zeta):
		return isSingletonClustering(G._this, zeta._this)
	@staticmethod
	def isOneClustering(Graph G, Partition zeta):
		return isOneClustering(G._this, zeta._this)
	@staticmethod
	def equalClustering(Partition zeta, Partition eta, Graph G):
		return equalClusterings(zeta._this, eta._this, G._this)

cdef extern from "cpp/graph/GraphTools.h" namespace "NetworKit::GraphTools":
	_Graph getCompactedGraph(_Graph G, unordered_map[node,node]) nogil except +
	unordered_map[node,node] getContinuousNodeIds(_Graph G) nogil except +
	unordered_map[node,node] getRandomContinuousNodeIds(_Graph G) nogil except +

cdef class GraphTools:
	@staticmethod
	def getCompactedGraph(Graph graph, nodeIdMap):
		"""
			Computes a graph with the same structure but with continuous node ids.
		"""
		cdef unordered_map[node,node] cNodeIdMap
		for key in nodeIdMap:
			cNodeIdMap[key] = nodeIdMap[key]
		return Graph().setThis(getCompactedGraph(graph._this,cNodeIdMap))

	@staticmethod
	def getContinuousNodeIds(Graph graph):
		"""
			Computes a map of node ids to continuous node ids.
		"""
		cdef unordered_map[node,node] cResult
		with nogil:
			cResult = getContinuousNodeIds(graph._this)
		result = dict()
		for elem in cResult:
			result[elem.first] = elem.second
		return result

	@staticmethod
	def getRandomContinuousNodeIds(Graph graph):
		"""
			Computes a map of node ids to continuous, randomly permutated node ids.
		"""
		cdef unordered_map[node,node] cResult
		with nogil:
			cResult = getRandomContinuousNodeIds(graph._this)
		result = dict()
		for elem in cResult:
			result[elem.first] = elem.second
		return result


cdef extern from "cpp/community/PartitionIntersection.h":
	cdef cppclass _PartitionIntersection "NetworKit::PartitionIntersection":
		_PartitionIntersection() except +
		_Partition calculate(_Partition zeta, _Partition eta) except +

cdef class PartitionIntersection:
	""" Class for calculating the intersection of two partitions, i.e. the clustering with the fewest clusters
	such that each cluster is a subset of a cluster in both partitions.
	"""
	cdef _PartitionIntersection _this
	def calculate(self, Partition zeta, Partition eta):
		"""  Calculate the intersection of two partitions `zeta` and `eta`

		Parameters
		----------
		zeta: Partition
			The first partition
		eta: Partition
			The second partition

		Returns
		-------
		Partition
			The intersection of zeta and eta
		"""
		return Partition().setThis(self._this.calculate(zeta._this, eta._this))

cdef extern from "cpp/community/Coverage.h":
	cdef cppclass _Coverage "NetworKit::Coverage":
		_Coverage() except +
		double getQuality(_Partition _zeta, _Graph _G) except +

cdef class Coverage:
	""" Coverage is the fraction of intra-community edges """
	cdef _Coverage _this

	def getQuality(self, Partition zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)


cdef extern from "cpp/community/EdgeCut.h":
	cdef cppclass _EdgeCut "NetworKit::EdgeCut":
		_EdgeCut() except +
		double getQuality(_Partition _zeta, _Graph _G) except +

cdef class EdgeCut:
	""" Edge cut is the total weight of inter-community edges"""
	cdef _EdgeCut _this

	def getQuality(self, Partition zeta, Graph G):
		return self._this.getQuality(zeta._this, G._this)


cdef extern from "cpp/community/Modularity.h":
	cdef cppclass _Modularity "NetworKit::Modularity":
		_Modularity() except +
		double getQuality(_Partition _zeta, _Graph _G) nogil except +


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
		cdef double ret
		with nogil:
			ret = self._this.getQuality(zeta._this, G._this)
		return ret

cdef extern from "cpp/community/HubDominance.h":
	cdef cppclass _HubDominance "NetworKit::HubDominance":
		_HubDominance() except +
		double getQuality(_Partition _zeta, _Graph _G) except +
		double getQuality(_Cover _zeta, _Graph _G) except +

cdef class HubDominance:
	"""
	A quality measure that measures the dominance of hubs in clusters. The hub dominance of a single
	cluster is defined as the maximum cluster-internal degree of a node in that cluster divided by
	the maximum cluster-internal degree, i.e. the number of nodes in the cluster minus one. The
	value for all clusters is defined as the average of all clusters.

	Strictly speaking this is not a quality measure as this is rather dependent on the type of the
	considered graph, for more information see
	Lancichinetti A, Kivelä M, Saramäki J, Fortunato S (2010)
	Characterizing the Community Structure of Complex Networks
	PLoS ONE 5(8): e11976. doi: 10.1371/journal.pone.0011976
	http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011976
	"""

	cdef _HubDominance _this

	def getQuality(self, PartitionCover zeta, Graph G):
		"""
		Calculates the dominance of hubs in the given Partition or Cover of the given
		Graph.

		Parameters
		----------
		zeta : Partition or Cover
			The Partition or Cover for which the hub dominance shall be calculated
		G : Graph
			The Graph to which zeta belongs

		Returns
		-------
		double
			The average hub dominance in the given Partition or Cover
		"""
		return self._this.getQuality(zeta._this, G._this)


cdef extern from "cpp/community/CommunityDetectionAlgorithm.h":
	cdef cppclass _CommunityDetectionAlgorithm "NetworKit::CommunityDetectionAlgorithm":
		_CommunityDetectionAlgorithm() # Workaround for Cython < 0.22
		_CommunityDetectionAlgorithm(const _Graph &_G)
		void run() nogil except +
		_Partition getPartition() except +
		string toString() except +


cdef class CommunityDetector:
	""" Abstract base class for static community detection algorithms """
	cdef _CommunityDetectionAlgorithm *_this
	cdef Graph _G

	def __init__(self, *args, **namedargs):
		if type(self) == CommunityDetector:
			raise RuntimeError("Error, you may not use CommunityDetector directly, use a sub-class instead")

	def __cinit__(self, *args, **namedargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL
		self._G = None # just to be sure the graph is deleted

	def run(self):
		"""
		Executes the community detection algorithm.

		Returns
		-------
		CommunityDetector:
			self
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		with nogil:
			self._this.run()
		return self

	def getPartition(self):
		"""  Returns a partition of the clustering.

		Returns
		-------
		Partition:
			A Partition of the clustering.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return Partition().setThis(self._this.getPartition())

	def toString(self):
		""" Get string representation.

		Returns
		-------
		string
			String representation of algorithm and parameters.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.toString().decode("utf-8")

cdef extern from "cpp/community/PLP.h":
	cdef cppclass _PLP "NetworKit::PLP"(_CommunityDetectionAlgorithm):
		_PLP(_Graph _G) except +
		_PLP(_Graph _G, count updateThreshold) except +
		_PLP(_Graph _G, _Partition baseClustering, count updateThreshold) except +
		_PLP(_Graph _G, _Partition baseClustering) except +
		count numberOfIterations() except +
		vector[count] getTiming() except +


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

	def __cinit__(self, Graph G not None, Partition baseClustering=None, updateThreshold=None):
		"""
		Constructor to the Parallel label propagation community detection algorithm.

		Parameters
		----------
		G : Graph
			The graph on which the algorithm has to run.
		baseClustering : Partition
			PLP needs a base clustering to start from; if none is given the algorithm will run on a singleton clustering.
		updateThreshold : integer
			number of nodes that have to be changed in each iteration so that a new iteration starts.
		"""
		self._G = G

		if updateThreshold is None and baseClustering is None:
			self._this = new _PLP(G._this)
		elif updateThreshold is None and baseClustering is not None:
			self._this = new _PLP(G._this, baseClustering._this)
		elif updateThreshold is not None and baseClustering is None:
			p = Partition(0)
			self._this = new _PLP(G._this, p._this, updateThreshold)
		else:
			self._this = new _PLP(G._this, baseClustering._this, updateThreshold)

	def numberOfIterations(self):
		""" Get number of iterations in last run.

		Returns
		-------
		count
			The number of iterations.
		"""
		return (<_PLP*>(self._this)).numberOfIterations()

	def getTiming(self):
		""" Get list of running times for each iteration.

		Returns
		-------
		count
			The list of running times in milliseconds.
		"""
		return (<_PLP*>(self._this)).getTiming()

cdef extern from "cpp/community/LPDegreeOrdered.h":
	cdef cppclass _LPDegreeOrdered "NetworKit::LPDegreeOrdered"(_CommunityDetectionAlgorithm):
		_LPDegreeOrdered(_Graph _G) except +
		count numberOfIterations()

cdef class LPDegreeOrdered(CommunityDetector):
	""" Label propagation-based community detection algorithm which processes nodes in increasing order of node degree.	"""

	def __cinit__(self, Graph G not None):
		self._G = G
		self._this = new _LPDegreeOrdered(G._this)

	def numberOfIterations(self):
		""" Get number of iterations in last run.

		Returns
		-------
		count
			Number of iterations.
		"""
		return (<_LPDegreeOrdered*>(self._this)).numberOfIterations()



cdef extern from "cpp/community/PLM.h":
	cdef cppclass _PLM "NetworKit::PLM"(_CommunityDetectionAlgorithm):
		_PLM(_Graph _G) except +
		_PLM(_Graph _G, bool refine, double gamma, string par, count maxIter, bool turbo) except +
		map[string, vector[count]] getTiming() except +

cdef extern from "cpp/community/PLM.h" namespace "NetworKit::PLM":
	pair[_Graph, vector[node]] PLM_coarsen "NetworKit::PLM::coarsen" (const _Graph& G, const _Partition& zeta) except +
	_Partition PLM_prolong "NetworKit::PLM::prolong"(const _Graph& Gcoarse, const _Partition& zetaCoarse, const _Graph& Gfine, vector[node] nodeToMetaNode) except +


cdef class PLM(CommunityDetector):
	""" Parallel Louvain Method - the Louvain method, optionally extended to
		a full multi-level algorithm with refinement

		Parameters
		----------
		G : Graph
			A graph.
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
		turbo : bool, optional
			faster but uses O(n) additional memory per thread
	"""

	def __cinit__(self, Graph G not None, refine=False, gamma=1.0, par="balanced", maxIter=32, turbo=False):
		self._G = G
		self._this = new _PLM(G._this, refine, gamma, stdstring(par), maxIter, turbo)

	def getTiming(self):
		"""  Get detailed time measurements.
		"""
		return (<_PLM*>(self._this)).getTiming()

	@staticmethod
	def coarsen(Graph G, Partition zeta, bool parallel = False):
		cdef pair[_Graph, vector[node]] result = move(PLM_coarsen(G._this, zeta._this))
		return (Graph().setThis(result.first), result.second)

	@staticmethod
	def prolong(Graph Gcoarse, Partition zetaCoarse, Graph Gfine, vector[node] nodeToMetaNode):
		return Partition().setThis(PLM_prolong(Gcoarse._this, zetaCoarse._this, Gfine._this, nodeToMetaNode))

cdef extern from "cpp/community/CutClustering.h":
	cdef cppclass _CutClustering "NetworKit::CutClustering"(_CommunityDetectionAlgorithm):
		_CutClustering(_Graph _G) except +
		_CutClustering(_Graph _G, edgeweight alpha) except +

cdef extern from "cpp/community/CutClustering.h" namespace "NetworKit::CutClustering":
	map[double, _Partition] CutClustering_getClusterHierarchy "NetworKit::CutClustering::getClusterHierarchy"(const _Graph& G) nogil except +


cdef class CutClustering(CommunityDetector):
	"""
	Cut clustering algorithm as defined in
	Flake, Gary William; Tarjan, Robert E.; Tsioutsiouliklis, Kostas. Graph Clustering and Minimum Cut Trees.
	Internet Mathematics 1 (2003), no. 4, 385--408.

	Parameters
	----------
	G : Graph
	alpha : double
		The parameter for the cut clustering algorithm
	"""
	def __cinit__(self, Graph G not None,  edgeweight alpha):
		self._G = G
		self._this = new _CutClustering(G._this, alpha)

	@staticmethod
	def getClusterHierarchy(Graph G not None):
		""" Get the complete hierarchy with all possible parameter values.

		Each reported parameter value is the lower bound for the range in which the corresponding clustering is calculated by the cut clustering algorithm.

		Warning: all reported parameter values are slightly too high in order to avoid wrong clusterings because of numerical inaccuracies.
		Furthermore the completeness of the hierarchy cannot be guaranteed because of these inaccuracies.
		This implementation hasn't been optimized for performance.

		Parameters
		----------
		G : Graph
			The graph.

		Returns
		-------
		dict
			A dictionary with the parameter values as keys and the corresponding Partition instances as values
		"""
		cdef map[double, _Partition] result
		# FIXME: this probably copies the whole hierarchy because of exception handling, using move might fix this
		with nogil:
			result = CutClustering_getClusterHierarchy(G._this)
		pyResult = {}
		# FIXME: this code copies the partitions a lot!
		for res in result:
			pyResult[res.first] = Partition().setThis(res.second)
		return pyResult

cdef class DissimilarityMeasure:
	""" Abstract base class for partition/community dissimilarity measures """
	# TODO: use conventional class design of parametrized constructor, run-method and getters
	pass


cdef extern from "cpp/community/NodeStructuralRandMeasure.h":
	cdef cppclass _NodeStructuralRandMeasure "NetworKit::NodeStructuralRandMeasure":
		_NodeStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class NodeStructuralRandMeasure(DissimilarityMeasure):
	""" The node-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering all pairs of nodes.
	"""
	cdef _NodeStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret


cdef extern from "cpp/community/GraphStructuralRandMeasure.h":
	cdef cppclass _GraphStructuralRandMeasure "NetworKit::GraphStructuralRandMeasure":
		_GraphStructuralRandMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class GraphStructuralRandMeasure(DissimilarityMeasure):
	""" The graph-structural Rand measure assigns a similarity value in [0,1]
		to two partitions of a graph, by considering connected pairs of nodes.
	"""
	cdef _GraphStructuralRandMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret


cdef extern from "cpp/community/JaccardMeasure.h":
	cdef cppclass _JaccardMeasure "NetworKit::JaccardMeasure":
		_JaccardMeasure() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class JaccardMeasure(DissimilarityMeasure):
	""" TODO:
	"""
	cdef _JaccardMeasure _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret

cdef extern from "cpp/community/NMIDistance.h":
	cdef cppclass _NMIDistance "NetworKit::NMIDistance":
		_NMIDistance() except +
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class NMIDistance(DissimilarityMeasure):
	""" The NMI distance assigns a similarity value in [0,1] to two partitions
		of a graph.
	"""
	cdef _NMIDistance _this

	def getDissimilarity(self, Graph G, Partition first, Partition second):
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret

cdef extern from "cpp/community/AdjustedRandMeasure.h":
	cdef cppclass _AdjustedRandMeasure "NetworKit::AdjustedRandMeasure":
		double getDissimilarity(_Graph G, _Partition first, _Partition second) nogil except +

cdef class AdjustedRandMeasure(DissimilarityMeasure):
	"""
	The adjusted rand dissimilarity measure as proposed by Huber and Arabie in "Comparing partitions" (http://link.springer.com/article/10.1007/BF01908075)
	"""
	cdef _AdjustedRandMeasure _this

	def getDissimilarity(self, Graph G not None, Partition first not None, Partition second not None):
		"""
		Get the adjust rand dissimilarity. Runs in O(n log(n)).

		Note that the dissimilarity can be larger than 1 if the partitions are more different than expected in the random model.

		Parameters
		----------
		G : Graph
			The graph on which the partitions shall be compared
		zeta : Partition
			The first partiton
		eta : Partition
			The second partition

		Returns
		-------
		double
			The adjusted rand dissimilarity
		"""
		cdef double ret
		with nogil:
			ret = self._this.getDissimilarity(G._this, first._this, second._this)
		return ret

cdef extern from "cpp/community/EPP.h":
	cdef cppclass _EPP "NetworKit::EPP"(_CommunityDetectionAlgorithm):
		_EPP(_Graph G)
		_Partition getCorePartition() except +
		vector[_Partition] getBasePartitions() except +

cdef class EPP(CommunityDetector):
	""" EPP - Ensemble Preprocessing community detection algorithm.
	Combines multiple base algorithms and a final algorithm. A consensus of the
	solutions of the base algorithms is formed and the graph is coarsened accordingly.
	Then the final algorithm operates on the coarse graph and determines a solution
	for the input graph.
	"""
	def __cinit__(self, Graph G not None):
		self._G = G
		self._this = new _EPP(G._this)

	def getCorePartition(self):
		"""  Returns the core partition the algorithm.

		Returns
		-------
		Partition:
			A Partition of the clustering.
		"""
		return Partition().setThis((<_EPP*>(self._this)).getCorePartition())

	def getBasePartitions(self):
		"""  Returns the base partitions of the algorithm.
		"""
		base = (<_EPP*>(self._this)).getBasePartitions()
		return [Partition().setThis(b) for b in base]

	cdef setThis(self, _EPP* other):
		del self._this # is this correct here?
		self._this = other
		return self

cdef extern from "cpp/community/EPPFactory.h" namespace "NetworKit::EPPFactory":
		#_EPP make(_Graph G, count ensembleSize, string baseAlgorithm, string finalAlgorithm)
		_EPP* makePtr(_Graph G, count ensembleSize, string baseAlgorithm, string finalAlgorithm)

cdef class EPPFactory:
	""" This class makes instaces of the EPP community detection algorithm """

	@staticmethod
	def make(Graph G not None, ensembleSize, baseAlgorithm="PLP", finalAlgorithm="PLM"):
		"""
		Returns an instance of an ensemble preprocessing (EPP).

		Parameters:
		-----------
		G : Graph
			The graph on which the ensemble is supposed to run.
		ensembleSize : integer
			The amount of baseAlgorithms to preprocess the communities.
		baseAlgorithm : CommunityDetectionAlgorithm
			String representation of the algorithm ("PLP","PLM") to preprocess the communities. ensembleSize instances will be created.
		finalAlgorithm  : CommunityDetectionAlgorithm
			String representation of the algorithm ("PLP" "PLM[R]") to finish the ensemble.

		Returns
		-------
		EPP
			The EPP instance.
		"""
		return EPP(G).setThis(makePtr(G._this, ensembleSize, stdstring(baseAlgorithm), stdstring(finalAlgorithm)))

# Module: flows

cdef extern from "cpp/flow/EdmondsKarp.h":
	cdef cppclass _EdmondsKarp "NetworKit::EdmondsKarp":
		_EdmondsKarp(const _Graph &graph, node source, node sink) except +
		void run() nogil except +
		edgeweight getMaxFlow() const
		vector[node] getSourceSet() except +
		edgeweight getFlow(node u, node v) except +
		edgeweight getFlow(edgeid eid) const
		vector[edgeweight] getFlowVector() except +

cdef class EdmondsKarp:
	"""
	The EdmondsKarp class implements the maximum flow algorithm by Edmonds and Karp.

	Parameters
	----------
	graph : Graph
		The graph
	source : node
		The source node for the flow calculation
	sink : node
		The sink node for the flow calculation
	"""
	cdef _EdmondsKarp* _this
	cdef Graph _graph

	def __cinit__(self, Graph graph not None, node source, node sink):
		self._graph = graph # store reference of graph for memory management, so the graph is not deallocated before this object
		self._this = new _EdmondsKarp(graph._this, source, sink)

	def __dealloc__(self):
		del self._this

	def run(self):
		"""
		Computes the maximum flow, executes the EdmondsKarp algorithm
		"""
		with nogil:
			self._this.run()
		return self

	def getMaxFlow(self):
		"""
		Returns the value of the maximum flow from source to sink.

		Returns
		-------
		edgeweight
			The maximum flow value
		"""
		return self._this.getMaxFlow()

	def getSourceSet(self):
		"""
		Returns the set of the nodes on the source side of the flow/minimum cut.

		Returns
		-------
		list
			The set of nodes that form the (smallest) source side of the flow/minimum cut.
		"""
		return self._this.getSourceSet()

	def getFlow(self, node u, node v = none):
		"""
		Get the flow value between two nodes u and v or an edge identified by the edge id u.
		Warning: The variant with two edge ids is linear in the degree of u.

		Parameters
		----------
		u : node or edgeid
			The first node incident to the edge or the edge id
		v : node
			The second node incident to the edge (optional if edge id is specified)

		Returns
		-------
		edgeweight
			The flow on the specified edge
		"""
		if v == none: # Assume that node and edge ids are the same type
			return self._this.getFlow(u)
		else:
			return self._this.getFlow(u, v)

	def getFlowVector(self):
		"""
		Return a copy of the flow values of all edges.

		Returns
		-------
		list
			The flow values of all edges indexed by edge id
		"""
		return self._this.getFlowVector()

# Module: properties

# this is an example for using static methods
cdef extern from "cpp/properties/GraphProperties.h" namespace "NetworKit::GraphProperties":
	# static methods live in the class namespace, so declare them here
	pair[count, count] minMaxDegree(_Graph _G) except +
	pair[pair[count, count], pair[count, count]] minMaxDegreeDirected(_Graph _G) except +
	double averageDegree(_Graph _G) except +
	vector[count] degreeDistribution(_Graph _G) except +
	vector[count] degreeSequence(_Graph _G) except +
	vector[double] localClusteringCoefficients(_Graph _G) except +
	double averageLocalClusteringCoefficient(_Graph _G) except +
	vector[double] localClusteringCoefficientPerDegree(_Graph _G) except +
	double degreeAssortativity(_Graph G, bool) except +
	double degreeAssortativityDirected(_Graph G, bool, bool) except +

	cdef cppclass _GraphProperties "NetworKit::GraphProperties":
		pass

cdef class GraphProperties:
	""" Collects various functions for basic graph properties """

	# TODO: this class should become obsolete with the new profiling module

	@staticmethod
	def minMaxDegree(Graph G not None):
		return minMaxDegree(G._this)

	@staticmethod
	def minMaxDegreeDirected(Graph G not None):
		return minMaxDegreeDirected(G._this)

	@staticmethod
	def averageDegree(Graph G not None):
		return averageDegree(G._this)

	@staticmethod
	def degreeDistribution(Graph G not None):
		return degreeDistribution(G._this)

	@staticmethod
	def degreeSequence(Graph G not None):
		return degreeSequence(G._this)

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
		return averageLocalClusteringCoefficient(G._this)

	@staticmethod
	def degreeAssortativity(Graph G, bool useWeights):
		""" Get degree assortativity of the undirected graph `G`.

		Parameters
		----------
		G : Graph
			The undirected graph
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
		return degreeAssortativity(G._this, useWeights)

	@staticmethod
	def degreeAssortativityDirected(Graph G, bool alpha, bool beta):
		""" Get degree assortativity of the directed graph `G`.

		Parameters
		----------
		G : Graph
			The directed graph
		alpha : bool
			When iterating over edges (u,v), if alpha is True, the out degree of u will be considered.
		beta : bool
			When iterating over edges (u,v), if beta is True, the out degree of v will be considered.

		Returns
		-------
		double
			Degree assortativity of the graph `G`.

		Notes
		-----
		Degree assortativity for directed graphs based on http://www.pnas.org/content/107/24/10815.full
		"""
		return degreeAssortativityDirected(G._this, alpha, beta)




cdef extern from "cpp/properties/ConnectedComponents.h":
	cdef cppclass _ConnectedComponents "NetworKit::ConnectedComponents":
		_ConnectedComponents(_Graph G) except +
		void run() nogil except +
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
	cdef Graph _G

	def __cinit__(self,  Graph G):
		self._G = G
		self._this = new _ConnectedComponents(G._this)

	def __dealloc__(self):
		del self._this

	def run(self):
		""" This method determines the connected components for the graph given in the constructor. """
		with nogil:
			self._this.run()
		return self

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


cdef extern from "cpp/properties/ParallelConnectedComponents.h":
	cdef cppclass _ParallelConnectedComponents "NetworKit::ParallelConnectedComponents":
		_ParallelConnectedComponents(_Graph G, bool coarsening) except +
		void run() nogil except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		_Partition getPartition() except +


cdef class ParallelConnectedComponents:
	""" Determines the connected components and associated values for
		an undirected graph.
	"""
	cdef _ParallelConnectedComponents* _this
	cdef Graph _G

	def __cinit__(self,  Graph G, coarsening=True	):
		self._G = G
		self._this = new _ParallelConnectedComponents(G._this, coarsening)

	def __dealloc__(self):
		del self._this

	def run(self):
		with nogil:
			self._this.run()
		return self

	def getPartition(self):
		return Partition().setThis(self._this.getPartition())

	def numberOfComponents(self):
		return self._this.numberOfComponents()

	def componentOfNode(self, v):
		return self._this.componentOfNode(v)


cdef extern from "cpp/properties/StronglyConnectedComponents.h":
	cdef cppclass _StronglyConnectedComponents "NetworKit::StronglyConnectedComponents":
		_StronglyConnectedComponents(_Graph G) except +
		void run() nogil except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		_Partition getPartition() except +


cdef class StronglyConnectedComponents:
	""" Determines the connected components and associated values for
		a directed graph.
	"""
	cdef _StronglyConnectedComponents* _this
	cdef Graph _G

	def __cinit__(self,  Graph G):
		self._G = G
		self._this = new _StronglyConnectedComponents(G._this)

	def __dealloc__(self):
		del self._this

	def run(self):
		with nogil:
			self._this.run()
		return self

	def getPartition(self):
		return Partition().setThis(self._this.getPartition())

	def numberOfComponents(self):
		return self._this.numberOfComponents()

	def componentOfNode(self, v):
		return self._this.componentOfNode(v)



cdef extern from "cpp/properties/ClusteringCoefficient.h" namespace "NetworKit::ClusteringCoefficient":
		double avgLocal(_Graph G) nogil except +
		double sequentialAvgLocal(_Graph G) nogil except +
		double approxAvgLocal(_Graph G, count trials) nogil except +
		double exactGlobal(_Graph G) nogil except +
		double approxGlobal(_Graph G, count trials) nogil except +

cdef class ClusteringCoefficient:
	@staticmethod
	def avgLocal(Graph G):
		"""
		DEPRECATED: Use centrality.LocalClusteringCoefficient and take average.

		This calculates the average local clustering coefficient of graph `G`.

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
		cdef double ret
		with nogil:
			ret = avgLocal(G._this)
		return ret

	@staticmethod
	def sequentialAvgLocal(Graph G):
		""" This calculates the average local clustering coefficient of graph `G` using inherently sequential triangle counting.
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
		cdef double ret
		with nogil:
			ret = sequentialAvgLocal(G._this)
		return ret

	@staticmethod
	def approxAvgLocal(Graph G, count trials):
		cdef double ret
		with nogil:
			ret = approxAvgLocal(G._this, trials)
		return ret

	@staticmethod
	def exactGlobal(Graph G):
		""" This calculates the global clustering coefficient. """
		cdef double ret
		with nogil:
			ret = exactGlobal(G._this)
		return ret

	@staticmethod
	def approxGlobal(Graph G, count trials):
		cdef double ret
		with nogil:
			ret = approxGlobal(G._this, trials)
		return ret


cdef extern from "cpp/properties/Diameter.h" namespace "NetworKit::Diameter":
	pair[count, count] estimatedDiameterRange(_Graph G, double error, pair[node,node] *proof) nogil except +
	count exactDiameter(_Graph G) nogil except +
	edgeweight estimatedVertexDiameter(_Graph G, count) nogil except +
	edgeweight estimatedVertexDiameterPedantic(_Graph G) nogil except +

cdef class Diameter:
	"""
	TODO: docstring
	"""

	@staticmethod
	def estimatedDiameterRange(Graph G, double error=0.1):
		""" Estimates a range for the diameter of @a G. Based on the algorithm suggested in
		C. Magnien, M. Latapy, M. Habib: Fast Computation of Empirically Tight Bounds for
		the Diameter of Massive Graphs. Journal of Experimental Algorithmics, Volume 13, Feb 2009.

		Returns
		-------
		pair
			Pair of lower and upper bound for diameter.
		"""
		cdef pair[count, count] ret
		with nogil:
			ret = estimatedDiameterRange(G._this, error, NULL)
		return ret

	@staticmethod
	def estimatedDiameterRangeWithProof(Graph G, double error=0.1):
		""" Estimates a range for the diameter of @a G like estimatedDiameterRange but also
		returns a pair of nodes that has the reported lower bound as distance if the graph is non-trivial.

		Returns
		-------
		tuple
			Tuple of two tuples, the first is the result and contains lower and upper bound, the second contains the two nodes whose distance is the reported lower bound
		"""
		cdef pair[node, node] proof
		cdef pair[count, count] result
		with nogil:
			result = estimatedDiameterRange(G._this, error, &proof)
		return (result, proof)

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
		cdef edgeweight diam
		with nogil:
			diam = exactDiameter(G._this)
		return diam

	@staticmethod
	def estimatedVertexDiameter(Graph G, count samples):
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
		cdef edgeweight diam
		with nogil:
			if samples == 0:
				diam = estimatedVertexDiameterPedantic(G._this)
			else:
				diam = estimatedVertexDiameter(G._this, samples)
		return diam

cdef extern from "cpp/properties/Eccentricity.h" namespace "NetworKit::Eccentricity":
	pair[node, count] getValue(_Graph G, node v) except +

cdef class Eccentricity:
	"""
	TODO: docstring
	"""

	@staticmethod
	def getValue(Graph G, v):
		return getValue(G._this, v)




cdef extern from "cpp/properties/EffectiveDiameter.h" namespace "NetworKit::EffectiveDiameter":
	double effectiveDiameter (_Graph G, double ratio, count k, count r) nogil except +
	double effectiveDiameterExact(_Graph G, double ratio) nogil except +
	map[count, double] hopPlot(_Graph G, count maxDistance, count k, count r) except +

cdef class EffectiveDiameter:

	@staticmethod
	def effectiveDiameter(Graph G, double ratio=0.9, count k=64, count r=7):
		""" Estimates the number of edges on average needed to reach 90% of all other nodes with a variaton of the ANF algorithm presented in the paper A Fast and Scalable Tool for Data Mining
			in Massive Graphs by Palmer, Gibbons and Faloutsos
		Parameters
		----------
		G : Graph
			The graph.
		ratio : double
			The percentage of nodes that shall be within stepwith
		k : count
			number of parallel approximations, bigger k -> longer runtime, more precise result
		r : count
			number of additional bits, important in tiny graphs
		Returns
		-------
		double
			the estimated effective diameter
		"""
		cdef double diam
		with nogil:
			diam = effectiveDiameter(G._this, ratio, k, r)
		return diam

	@staticmethod
	def effectiveDiameterExact(Graph G, double ratio=0.9):
		""" Calculates the number of edges on average needed to reach 90% of all other nodes
		Parameters
		----------
		G : Graph
			The graph.
		ratio : double
			The percentage of nodes that shall be within stepwith
		Returns
		-------
		double
			the effective diameter
		"""
		cdef double diam
		with nogil:
			diam = effectiveDiameterExact(G._this, ratio)
		return diam

	@staticmethod
	def hopPlot(Graph G, maxDistance=0, k=64, r=7):
		""" Calculates the number of connected nodes for each distance between 0 and the diameter of the graph
		Parameters
		----------
		G : Graph
			The graph.
		maxDistance : double
			maximum distance between considered nodes
			set to 0 or negative to get the hop-plot for the entire graph so that each node can reach each other node
		k : count
			number of parallel approximations, bigger k -> longer runtime, more precise result
		r : count
			number of additional bits, important in tiny graphs
		Returns
		-------
		map
			number of connected nodes for each distance
		"""
		return hopPlot(G._this, maxDistance, k, r)


# Module: centrality

cdef extern from "cpp/centrality/Centrality.h":
	cdef cppclass _Centrality "NetworKit::Centrality":
		_Centrality(_Graph, bool, bool) except +
		void run() nogil except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +
		double maximum() except +


cdef class Centrality:
	""" Abstract base class for centrality measures"""

	cdef _Centrality* _this
	cdef Graph _G

	def __init__(self, *args, **kwargs):
		if type(self) == Centrality:
			raise RuntimeError("Error, you may not use Centrality directly, use a sub-class instead")

	def __cinit__(self, *args, **kwargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL
		self._G = None # just to be sure the graph is deleted

	def run(self):
		"""
		Executes the centrality algorithm.

		Returns
		-------
		Centrality:
			self
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		with nogil:
			self._this.run()
		return self

	def scores(self):
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.scores()

	def score(self, v):
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.score(v)

	def ranking(self):
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.ranking()

	def maximum(self):
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.maximum()


cdef extern from "cpp/centrality/DegreeCentrality.h":
	cdef cppclass _DegreeCentrality "NetworKit::DegreeCentrality" (_Centrality):
		_DegreeCentrality(_Graph, bool normalized) except +

cdef class DegreeCentrality(Centrality):
	""" Node centrality index which ranks nodes by their degree.
 	Optional normalization by maximum degree.

 	DegreeCentrality(G, normalized=False)

 	Constructs the DegreeCentrality class for the given Graph `G`. If the scores should be normalized,
 	then set `normalized` to True.

 	Parameters
 	----------
 	G : Graph
 		The graph.
 	normalized : bool, optional
 		Normalize centrality values in the interval [0,1].
	"""

	def __cinit__(self, Graph G, bool normalized=False):
		self._G = G
		self._this = new _DegreeCentrality(G._this, normalized)



cdef extern from "cpp/centrality/Betweenness.h":
	cdef cppclass _Betweenness "NetworKit::Betweenness" (_Centrality):
		_Betweenness(_Graph, bool, bool) except +
		vector[double] edgeScores() except +

cdef class Betweenness(Centrality):
	"""
		Betweenness(G, normalized=False, computeEdgeCentrality=False)

		Constructs the Betweenness class for the given Graph `G`. If the betweenness scores should be normalized,
  		then set `normalized` to True.

	 	Parameters
	 	----------
	 	G : Graph
	 		The graph.
	 	normalized : bool, optional
	 		Set this parameter to True if scores should be normalized in the interval [0,1].
		computeEdgeCentrality: bool, optional
			Set this to true if edge betweenness scores should be computed as well.
	"""

	def __cinit__(self, Graph G, normalized=False, computeEdgeCentrality=False):
		self._G = G
		self._this = new _Betweenness(G._this, normalized, computeEdgeCentrality)


	def edgeScores(self):
		""" Get a vector containing the betweenness score for each edge in the graph.

		Returns
		-------
		vector
			The betweenness scores calculated by run().
		"""
		return (<_Betweenness*>(self._this)).edgeScores()


cdef extern from "cpp/centrality/Closeness.h":
	cdef cppclass _Closeness "NetworKit::Closeness" (_Centrality):
		_Closeness(_Graph, bool) except +

cdef class Closeness(Centrality):
	"""
		Closeness(G, normalized=False)

		Constructs the Closeness class for the given Graph `G`. If the Closeness scores should be normalized,
  		then set `normalized` to True.

	 	Parameters
	 	----------
	 	G : Graph
	 		The graph.
	 	normalized : bool, optional
	 		Set this parameter to True if scores should be normalized in the interval [0,1]. Normalization only for unweighted networks.
	"""

	def __cinit__(self, Graph G, normalized=False):
		self._G = G
		self._this = new _Closeness(G._this, normalized)


cdef extern from "cpp/centrality/KPathCentrality.h":
	cdef cppclass _KPathCentrality "NetworKit::KPathCentrality" (_Centrality):
		_KPathCentrality(_Graph, double, count) except +

cdef class KPathCentrality(Centrality):
	"""
		KPathCentrality(G, alpha=0.2, k=0)

		Constructs the K-Path Centrality class for the given Graph `G`.

	 	Parameters
	 	----------
	 	G : Graph
	 		The graph.
	 	alpha : double, in interval [-0.5, 0.5]
			tradeoff between runtime and precision
			-0.5: maximum precision, maximum runtime
	 		 0.5: lowest precision, lowest runtime
	"""

	def __cinit__(self, Graph G, alpha=0.2, k=0):
		self._G = G
		self._this = new _KPathCentrality(G._this, alpha, k)


cdef extern from "cpp/centrality/KatzCentrality.h":
	cdef cppclass _KatzCentrality "NetworKit::KatzCentrality" (_Centrality):
		_KatzCentrality(_Graph, double, count) except +

cdef class KatzCentrality(Centrality):
	"""
		KatzCentrality(G, alpha=5e-4, beta=0.1, tol=1e-8)

		Constructs a KatzCentrality object for the given Graph `G`

	 	Parameters
	 	----------
	 	G : Graph
	 		The graph.
	 	alpha : double
			Damping of the matrix vector product result
		beta : double
			Constant value added to the centrality of each vertex
		tol : double
			The tolerance for convergence.
	"""

	def __cinit__(self, Graph G, alpha=0.2, k=0):
		self._G = G
		self._this = new _KatzCentrality(G._this, alpha, k)




cdef extern from "cpp/centrality/ApproxBetweenness.h":
	cdef cppclass _ApproxBetweenness "NetworKit::ApproxBetweenness" (_Centrality):
		_ApproxBetweenness(_Graph, double, double, count) except +
		count numberOfSamples() except +

cdef class ApproxBetweenness(Centrality):
	""" Approximation of betweenness centrality according to algorithm described in
 	Matteo Riondato and Evgenios M. Kornaropoulos: Fast Approximation of Betweenness Centrality through Sampling

 	ApproxBetweenness(G, epsilon=0.01, delta=0.1)

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

	def __cinit__(self, Graph G, epsilon=0.01, delta=0.1, diameterSamples=0):
		self._G = G
		self._this = new _ApproxBetweenness(G._this, epsilon, delta, diameterSamples)

	def numberOfSamples(self):
		return (<_ApproxBetweenness*>(self._this)).numberOfSamples()



cdef extern from "cpp/centrality/ApproxBetweenness2.h":
	cdef cppclass _ApproxBetweenness2 "NetworKit::ApproxBetweenness2" (_Centrality):
		_ApproxBetweenness2(_Graph, count, bool) except +


cdef class ApproxBetweenness2(Centrality):
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

	def __cinit__(self, Graph G, nSamples, normalized=False):
		self._G = G
		self._this = new _ApproxBetweenness2(G._this, nSamples, normalized)


cdef extern from "cpp/centrality/PageRank.h":
	cdef cppclass _PageRank "NetworKit::PageRank" (_Centrality):
		_PageRank(_Graph, double damp, double tol) except +

cdef class PageRank(Centrality):
	"""	Compute PageRank as node centrality measure.

	PageRank(G, damp=0.85, tol=1e-9)

	Parameters
	----------
	G : Graph
		Graph to be processed.
	damp : double
		Damping factor of the PageRank algorithm.
	tol : double, optional
		Error tolerance for PageRank iteration.
	"""

	def __cinit__(self, Graph G, double damp=0.85, double tol=1e-9):
		self._G = G
		self._this = new _PageRank(G._this, damp, tol)



cdef extern from "cpp/centrality/EigenvectorCentrality.h":
	cdef cppclass _EigenvectorCentrality "NetworKit::EigenvectorCentrality" (_Centrality):
		_EigenvectorCentrality(_Graph, double tol) except +

cdef class EigenvectorCentrality(Centrality):
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

	def __cinit__(self, Graph G, double tol=1e-9):
		self._G = G
		self._this = new _EigenvectorCentrality(G._this, tol)


cdef extern from "cpp/centrality/CoreDecomposition.h":
	cdef cppclass _CoreDecomposition "NetworKit::CoreDecomposition" (_Centrality):
		_CoreDecomposition(_Graph)
		vector[set[node]] cores() except +
		vector[set[node]] shells() except +
		index maxCoreNumber() except +

cdef class CoreDecomposition(Centrality):
	""" Computes k-core decomposition of a graph.

	CoreDecomposition(G)

	Create CoreDecomposition class for graph `G`.

	Parameters
	----------
	G : Graph
		The graph.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _CoreDecomposition(G._this)

	def maxCoreNumber(self):
		""" Get maximum core number.

		Returns
		-------
		index
			The maximum core number.
		"""
		return (<_CoreDecomposition*>(self._this)).maxCoreNumber()

	def cores(self):
		""" Get the k-cores as sets of nodes, indexed by k.

		Returns
		-------
		vector
			The k-cores as sets of nodes, indexed by k.
		"""
		return (<_CoreDecomposition*>(self._this)).cores()

	def shells(self):
		""" Get the k-shells as sets of nodes, indexed by k.

		Returns
		-------
		vector
			The k-shells as sets of nodes, indexed by k.
		"""
		return (<_CoreDecomposition*>(self._this)).shells()


cdef extern from "cpp/centrality/LocalClusteringCoefficient.h":
	cdef cppclass _LocalClusteringCoefficient "NetworKit::LocalClusteringCoefficient" (_Centrality):
		_LocalClusteringCoefficient(_Graph) except +

cdef class LocalClusteringCoefficient(Centrality):
	"""
		LocalClusteringCoefficient(G, normalized=False, computeEdgeCentrality=False)

		Constructs the LocalClusteringCoefficient class for the given Graph `G`. If the local clustering coefficient values should be normalized,
  		then set `normalized` to True.

	 	Parameters
	 	----------
	 	G : Graph
	 		The graph.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _LocalClusteringCoefficient(G._this)


cdef extern from "cpp/centrality/DynApproxBetweenness.h":
	cdef cppclass _DynApproxBetweenness "NetworKit::DynApproxBetweenness":
		_DynApproxBetweenness(_Graph, double, double, bool) except +
		void run() nogil except +
		void update(vector[_GraphEvent]) except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +
		count getNumberOfSamples() except +

cdef class DynApproxBetweenness:
	""" New dynamic algorithm for the approximation of betweenness centrality with
	a guaranteed error

	DynApproxBetweenness(G, epsiolon=0.01, delta=0.1, [storePredecessors])

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
	storePredecessors : bool
		store lists of predecessors?
	"""
	cdef _DynApproxBetweenness* _this
	cdef Graph _G

	def __cinit__(self, Graph G, epsilon=0.01, delta=0.1, storePredecessors = True):
		self._G = G
		self._this = new _DynApproxBetweenness(G._this, epsilon, delta, storePredecessors)

	# this is necessary so that the C++ object gets properly garbage collected
	def __dealloc__(self):
		del self._this

	def run(self):
		with nogil:
			self._this.run()
		return self

	def update(self, batch):
		""" Updates the betweenness centralities after the batch `batch` of edge insertions.

		Parameters
		----------
		batch : list of GraphEvent.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.update(_batch)

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

	def getNumberOfSamples(self):
		"""
		Get number of path samples used in last calculation.
		"""
		return self._this.getNumberOfSamples()


# Module: distmeasures

cdef extern from "cpp/distmeasures/AlgebraicDistance.h":
	cdef cppclass _AlgebraicDistance "NetworKit::AlgebraicDistance":
		_AlgebraicDistance(const _Graph& G, count numberSystems, count numberIterations, double omega, index norm) except +
		void preprocess() nogil except +
		double distance(node u, node v) except +
		vector[vector[double]] getLoadsOnNodes()

cdef class AlgebraicDistance:
	"""
	Algebraic distance assigns a distance value to pairs of nodes
	according to their structural closeness in the graph.

	Parameters
	----------
	G : Graph
		The graph.
	numberSystems : count
		Number of vectors/systems used for algebraic iteration.
	numberIterations : count
		Number of iterations in each system.
	omega : double, optional
		Overrelaxation parameter, default: 0.5.
	norm : index, optional
		The norm factor of the extended algebraic distance. Maximum norm is realized by setting the norm to 0. Default: 2.
	"""
	cdef _AlgebraicDistance* _this
	cdef Graph _G

	def __cinit__(self, Graph G, count numberSystems, count numberIterations, double omega = 0.5, index norm = 2):
		self._G = G
		self._this = new _AlgebraicDistance(G._this, numberSystems, numberIterations, omega, norm)

	def __dealloc__(self):
		del self._this

	def preprocess(self):
		"""
		Starting with random initialization, compute for all numberSystems
		"diffusion" systems the situation after numberIterations iterations
		of overrelaxation with overrelaxation parameter omega.

		REQ: Needs to be called before algdist delivers meaningful results!
		"""
		with nogil:
			self._this.preprocess()


	def distance(self, node u, node v):
		"""
		Returns the extended algebraic distance between node u and node v in the norm specified in
		the constructor.

		Parameters
		----------
		u : node
			The first node
		v : node
			The second node

		Returns
		-------
		Extended algebraic distance between the two nodes.
		"""
		return self._this.distance(u, v)


	def getLoadsOnNodes(self):
		"""
		Returns a list, indexed by node id, of the load values.
		"""
		return self._this.getLoadsOnNodes()

# Module: dynamic

cdef extern from "cpp/dynamics/GraphEvent.h":
	enum _GraphEventType "NetworKit::GraphEvent::Type":
		NODE_ADDITION,
		NODE_REMOVAL,
		EDGE_ADDITION,
		EDGE_REMOVAL,
		EDGE_WEIGHT_UPDATE,
		TIME_STEP

cdef extern from "cpp/dynamics/GraphEvent.h":
	cdef cppclass _GraphEvent "NetworKit::GraphEvent":
		node u, v
		edgeweight w
		_GraphEventType type
		_GraphEvent() except +
		_GraphEvent(_GraphEventType type, node u, node v, edgeweight w) except +
		string toString() except +

cdef class GraphEvent:
	cdef _GraphEvent _this

	NODE_ADDITION = 0
	NODE_REMOVAL = 1
	EDGE_ADDITION = 2
	EDGE_REMOVAL = 3
	EDGE_WEIGHT_UPDATE = 4
	TIME_STEP = 5

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

	def __cinit__(self, _GraphEventType type, node u, node v, edgeweight w):
		self._this = _GraphEvent(type, u, v, w)

	def toString(self):
		return self._this.toString().decode("utf-8")

	def __repr__(self):
		return self.toString()


cdef extern from "cpp/dynamics/DGSStreamParser.h":
	cdef cppclass _DGSStreamParser "NetworKit::DGSStreamParser":
		_DGSStreamParser(string path, bool mapped, node baseIndex) except +
		vector[_GraphEvent] getStream() except +

cdef class DGSStreamParser:
	cdef _DGSStreamParser* _this

	def __cinit__(self, path, mapped=True, baseIndex=0):
		self._this = new _DGSStreamParser(stdstring(path), mapped, baseIndex)

	def __dealloc__(self):
		del self._this

	def getStream(self):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.getStream()]


cdef extern from "cpp/dynamics/DGSWriter.h":
	cdef cppclass _DGSWriter "NetworKit::DGSWriter":
		void write(vector[_GraphEvent] stream, string path) except +


cdef class DGSWriter:
	cdef _DGSWriter* _this

	def __cinit__(self):
		self._this = new _DGSWriter()

	def __dealloc__(self):
		del self._this

	def write(self, stream, path):
		cdef vector[_GraphEvent] _stream
		for ev in stream:
			_stream.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.write(_stream, stdstring(path))


# cdef extern from "cpp/dcd2/DynamicCommunityDetection.h":
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



cdef extern from "cpp/generators/DynamicPathGenerator.h":
	cdef cppclass _DynamicPathGenerator "NetworKit::DynamicPathGenerator":
		_DynamicPathGenerator() except +
		vector[_GraphEvent] generate(count nSteps) except +


cdef class DynamicPathGenerator:
	""" Example dynamic graph generator: Generates a dynamically growing path. """
	cdef _DynamicPathGenerator* _this

	def __cinit__(self):
		self._this = new _DynamicPathGenerator()

	def __dealloc__(self):
		del self._this

	def generate(self, nSteps):
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]


cdef extern from "cpp/generators/DynamicDorogovtsevMendesGenerator.h":
	cdef cppclass _DynamicDorogovtsevMendesGenerator "NetworKit::DynamicDorogovtsevMendesGenerator":
		_DynamicDorogovtsevMendesGenerator() except +
		vector[_GraphEvent] generate(count nSteps) except +


cdef class DynamicDorogovtsevMendesGenerator:
	""" Generates a graph according to the Dorogovtsev-Mendes model.

 	DynamicDorogovtsevMendesGenerator()

 	Constructs the generator class.
	"""
	cdef _DynamicDorogovtsevMendesGenerator* _this

	def __cinit__(self):
		self._this = new _DynamicDorogovtsevMendesGenerator()

	def __dealloc__(self):
		del self._this

	def generate(self, nSteps):
		""" Generate event stream.

		Parameters
		----------
		nSteps : count
			Number of time steps in the event stream.
		"""
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]



cdef extern from "cpp/generators/DynamicPubWebGenerator.h":
	cdef cppclass _DynamicPubWebGenerator "NetworKit::DynamicPubWebGenerator":
		_DynamicPubWebGenerator(count numNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors) except +
		vector[_GraphEvent] generate(count nSteps) except +
		_Graph getGraph() except +


cdef class DynamicPubWebGenerator:
	cdef _DynamicPubWebGenerator* _this

	def __cinit__(self, numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors):
		self._this = new _DynamicPubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	def __dealloc__(self):
		del self._this

	def generate(self, nSteps):
		""" Generate event stream.

		Parameters
		----------
		nSteps : count
			Number of time steps in the event stream.
		"""
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]

	def getGraph(self):
		return Graph().setThis(self._this.getGraph())





cdef extern from "cpp/generators/DynamicForestFireGenerator.h":
	cdef cppclass _DynamicForestFireGenerator "NetworKit::DynamicForestFireGenerator":
		_DynamicForestFireGenerator(double p, bool directed, double r) except +
		vector[_GraphEvent] generate(count nSteps) except +
		_Graph getGraph() except +


cdef class DynamicForestFireGenerator:
	""" Generates a graph according to the forest fire model.
	 The forest fire generative model produces dynamic graphs with the following properties:
     heavy tailed degree distribution
     communities
     densification power law
     shrinking diameter

    see Leskovec, Kleinberg, Faloutsos: Graphs over Tim: Densification Laws,
    Shringking Diameters and Possible Explanations

 	DynamicForestFireGenerator(double p, bool directed, double r = 1.0)

 	Constructs the generator class.

 	Parameters
 	----------
 	p : forward burning probability.
 	directed : decides whether the resulting graph should be directed
 	r : optional, backward burning probability
	"""
	cdef _DynamicForestFireGenerator* _this

	def __cinit__(self, p, directed, r = 1.0):
		self._this = new _DynamicForestFireGenerator(p, directed, r)

	def __dealloc__(self):
		del self._this

	def generate(self, nSteps):
		""" Generate event stream.

		Parameters
		----------
		nSteps : count
			Number of time steps in the event stream.
		"""
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]




cdef extern from "cpp/dynamics/GraphUpdater.h":
	cdef cppclass _GraphUpdater "NetworKit::GraphUpdater":
		_GraphUpdater(_Graph G) except +
		void update(vector[_GraphEvent] stream) nogil except +
		vector[pair[count, count]] getSizeTimeline() except +

cdef class GraphUpdater:
	""" Updates a graph according to a stream of graph events.

	Parameters
	----------
	G : Graph
	 	initial graph
	"""
	cdef _GraphUpdater* _this
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _GraphUpdater(G._this)

	def __dealloc__(self):
		del self._this

	def update(self, stream):
		cdef vector[_GraphEvent] _stream
		for ev in stream:
			_stream.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		with nogil:
			self._this.update(_stream)


# Module: coarsening

cdef extern from "cpp/coarsening/ParallelPartitionCoarsening.h":
	cdef cppclass _ParallelPartitionCoarsening "NetworKit::ParallelPartitionCoarsening":
		_ParallelPartitionCoarsening() except +
		pair[_Graph, vector[node]] run(_Graph, _Partition) except +


cdef class ParallelPartitionCoarsening:
	cdef _ParallelPartitionCoarsening* _this

	def __cinit__(self):
		self._this = new _ParallelPartitionCoarsening()

	def __dealloc__(self):
		del self._this

	def run(self, Graph G not None, Partition zeta not None):
		result = self._this.run(G._this, zeta._this)
		return (Graph(0).setThis(result.first), result.second)

# Module: scd

cdef extern from "cpp/scd/PageRankNibble.h":
	cdef cppclass _PageRankNibble "NetworKit::PageRankNibble":
		_PageRankNibble(_Graph G, double epsilon, double alpha) except +
		map[node, set[node]] run(set[unsigned int] seeds) except +

cdef class PageRankNibble:
	"""
	Produces a cut around a given seed node using the PageRank-Nibble algorithm.
	see Andersen, Chung, Lang: Local Graph Partitioning using PageRank Vectors

	Parameters:
	-----------
	G : graph in which the cut is to be produced, must be unweighted.
	epsilon : the max probability in the residual vector for each node.
	alpha : the random walk loop probability.
	"""
	cdef _PageRankNibble *_this
	cdef Graph _G

	def __cinit__(self, Graph G, double epsilon, double alpha):
		self._G = G
		self._this = new _PageRankNibble(G._this, epsilon, alpha)

	def run(self, set[unsigned int] seeds):
		"""
		Produces a cut around a given seed node.

		Parameters:
		-----------
		seeds : the seed node ids.
		"""
		return self._this.run(seeds)

cdef extern from "cpp/scd/GCE.h":
	cdef cppclass _GCE "NetworKit::GCE":
		_GCE(_Graph G, string quality) except +
		map[node, set[node]] run(set[unsigned int] seeds) except +

cdef class GCE:
	"""
	Produces a cut around a given seed node using the GCE algorithm.

	Parameters:
	-----------
	G : graph in which the cut is to be produced, must be unweighted.
	"""
	cdef _GCE *_this
	cdef Graph _G

	def __cinit__(self, Graph G, quality):
		self._G = G
		self._this = new _GCE(G._this, stdstring(quality))

	def run(self, set[unsigned int] seeds):
		"""
		Produces a cut around a given seed node.

		Parameters:
		-----------
		seeds : the seed node ids.
		"""
		return self._this.run(seeds)

# Module: clique

cdef extern from "cpp/clique/MaxClique.h":
	cdef cppclass _MaxClique "NetworKit::MaxClique":
		_MaxClique(_Graph G, count lb) except +
		void run() nogil except +
		count getMaxCliqueSize() except +

cdef class MaxClique:
	"""
	Exact algorithm for computing the size of the largest clique in a graph.
	Worst-case running time is exponential, but in practice the algorithm is fairly fast.
	Reference: Pattabiraman et al., http://arxiv.org/pdf/1411.7460.pdf

	Parameters:
	-----------
	G : graph in which the cut is to be produced, must be unweighted.
	lb : the lower bound of the size of the maximum clique.
	"""
	cdef _MaxClique* _this
	cdef Graph _G

	def __cinit__(self, Graph G not None, lb=0):
		self._G = G
		self._this = new _MaxClique(G._this, lb)


	def __dealloc__(self):
		del self._this

	def run(self):
		"""
		Actual maximum clique algorithm. Determines largest clique each vertex
	 	is contained in and returns size of largest. Pruning steps keep running time
	 	acceptable in practice.
	 	"""
		cdef count size
		with nogil:
			self._this.run()

	def getMaxCliqueSize(self):
		"""
		Returns the size of the biggest clique
		"""
		return self._this.getMaxCliqueSize()
