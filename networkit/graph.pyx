# distutils: language=c++

from cython.operator import dereference, preincrement

from .base import Algorithm
from .helpers import stdstring, pystring
from .traversal import Traversal
from . import graphio
import os

cdef class Graph:

	""" 
	Graph(n=0, weighted=False, directed=False, edgesIndexed=False)

	An undirected graph (with optional weights) and parallel iterator methods.

	Create a graph of `n` nodes. The graph has assignable edge weights if `weighted` is set to True.
	If `weighted` is set to False each edge has edge weight 1.0 and any other weight assignment will
	be ignored.

	Parameters
	----------
	n : int, optional
		Number of nodes.
	weighted : bool, optional
		If set to True, the graph can have edge weights other than 1.0. Default: False
	directed : bool, optional
		If set to True, the graph will be directed. Default: False
	edgesIndexed : bool, optional
		If set to True, the graph's edges will be indexed. Default: False
	"""

	def __cinit__(self, n=0, bool_t weighted=False, bool_t directed=False, bool_t edgesIndexed=False):
		if isinstance(n, Graph):
			self._this = move(_Graph((<Graph>n)._this, weighted, directed, edgesIndexed))
		else:
			self._this = move(_Graph(<count>n, weighted, directed, edgesIndexed))

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
		return "NetworKit.Graph(n={0}, m={1})".format(self.numberOfNodes(), self.numberOfEdges())
	
	def __getstate__(self):
		return graphio.NetworkitBinaryWriter(graphio.Format.NetworkitBinary, chunks = 32, weightsType = 5).writeToBuffer(self)
	
	def __setstate__(self, state):
		newG = graphio.NetworkitBinaryReader().readFromBuffer(state)
		self._this = move(_Graph((<Graph>newG)._this, <bool_t>(newG.isWeighted()), <bool_t>(newG.isDirected()), <bool_t>(newG.hasEdgeIds())))

	def indexEdges(self, bool_t force = False):
		"""
		indexEdges(force = False)

		Assign integer ids to edges.

		Parameters
		----------
		force : bool, optional
			Force re-indexing of edges.
		"""
		self._this.indexEdges(force)

	def hasEdgeIds(self):
		"""
		hasEdgeIds()

		Returns true if edges have been indexed

		Returns
		-------
		bool
			If edges have been indexed
		"""
		return self._this.hasEdgeIds()

	def edgeId(self, node u, node v):
		"""
		edgeId(u, v)

		Returns
		-------
		int
			Id of the edge.
		"""
		return self._this.edgeId(u, v)

	def numberOfNodes(self):
		"""
		numberOfNodes()

		Get the number of nodes in the graph.

		Returns
		-------
		int
			The number of nodes.
		"""
		return self._this.numberOfNodes()

	def numberOfEdges(self):
		"""
		numberOfEdges()

		Get the number of edges in the graph.

		Returns
	 	-------
		int
			The number of edges.
		"""
		return self._this.numberOfEdges()

	def upperNodeIdBound(self):
		"""
		upperNodeIdBound()

		Get an upper bound for the node ids in the graph.

		Returns
		-------
		int
			An upper bound for the node ids in the graph.
		"""
		return self._this.upperNodeIdBound()

	def upperEdgeIdBound(self):
		"""
		upperEdgeIdBound()

		Get an upper bound for the edge ids in the graph.

		Returns
		-------
		int
			An upper bound for the edge ids in the graph.
		"""
		return self._this.upperEdgeIdBound()

	def degree(self, u):
		"""
		degree(u)

		Get the number of neighbors of `v`.

		Parameters
		----------
		v : int
			The input Node.

		Returns
		-------
		int
			The number of neighbors.
		"""
		return self._this.degree(u)

	def degreeIn(self, u):
		"""
		degreeIn(u)

		Get the number of in-neighbors of `v`.

		Parameters
		----------
		v : int
			The input Node.

		Returns
		-------
		int
			The number of in-neighbors.
		"""		
		return self._this.degreeIn(u)

	def degreeOut(self, u):
		"""
		degreeOut(u)

		Get the number of out-neighbors of `v`.

		Parameters
		----------
		v : int
			The Input Node.i
		Returns
		-------
		int
			The number of out-neighbors.
		"""
		return self._this.degreeOut(u)

	def weightedDegree(self, u, countSelfLoopsTwice=False):
		"""
		weightedDegree(u, countSelfLoopsTwice=False)

		Returns the weighted out-degree of u.

		For directed graphs this is the sum of weights of all outgoing edges of u.

		Parameters
		----------
		u : int
			The input Node.
		countSelfLoopsTwice : bool, optional
			If set to True, self-loops will be counted twice. Default: False

		Returns
		-------
		float
			The weighted out-degree of u.
		"""
		return self._this.weightedDegree(u, countSelfLoopsTwice)

	def weightedDegreeIn(self, u, countSelfLoopsTwice=False):
		"""
		weightedDegreeIn(u, countSelfLoopsTwice=False)

		Returns the weighted in-degree of u.

		For directed graphs this is the sum of weights of all ingoing edges of u.

		Parameters
		----------
		u : int
			The input node.
		countSelfLoopsTwice : bool, optional
			If set to True, self-loops will be counted twice. Default: False

		Returns
		-------
		float
			The weighted in-degree of u.
		"""
		return self._this.weightedDegreeIn(u, countSelfLoopsTwice)

	def isIsolated(self, u):
		"""
		isIsolated(u)

		If the node `u` is isolated.

		Parameters
		----------
		u : int
			The input node.

		Returns
		-------
		bool
			Indicates whether the node is isolated.
		"""
		return self._this.isIsolated(u)

	def addNode(self):
		""" 
		addNode()
		
		Add a new node to the graph and return it.

		Returns
		-------
		int
			The new node.
		"""
		return self._this.addNode()

	def addNodes(self, numberOfNewNodes):
		""" 
		addNodes(numberOfNewNodes)
		
		Add numberOfNewNodes many new nodes to the graph and return
		the id of the last node added.

		Parameters
		----------
		numberOfNewNodes : int
			Number of nodes to be added.

		Returns
		-------
		int
			The id of the last node added.
		"""
		assert(numberOfNewNodes >= 0)
		return self._this.addNodes(numberOfNewNodes)

	def removeNode(self, u):
		"""
		removeNode(u)
		
		Remove a node `u` and all incident edges from the graph.

		Incoming as well as outgoing edges will be removed.

		Parameters
		----------
		u : int
			Id of node to be removed.
		"""
		self._this.removeNode(u)

	def restoreNode(self, u):
		""" 
		restoreNode(u)

		Restores a previously deleted node `u` with its previous id in the graph.

		Parameters
		----------
		u : int
			The input node.
		"""
		self._this.restoreNode(u)

	def hasNode(self, u):
		""" 
		hasNode(u)
		
		Checks if the Graph has the node `u`, i.e. if `u` hasn't been deleted and is in the range of valid ids.

		Parameters
		----------
		u : int
			Id of node queried.

		Returns
		-------
		bool
			Indicates whether node `u` is part of the graph.
		"""
		return self._this.hasNode(u)

	def addEdge(self, u, v, w=1.0, addMissing = False, checkMultiEdge = False):
		""" 
		addEdge(u, v, w=1.0, addMissing=False, checkMultiEdge=False)
		
		Insert an undirected edge between the nodes `u` and `v`. If the graph is weighted you can optionally set a weight for this edge. 
		The default weight is 1.0. If one or both end-points do not exists and addMissing is set, they are silently added.
		
		Note
		----
		By default it is not checked whether this edge already exists, thus it is possible to create multi-edges. Multi-edges are not supported and will NOT be
		handled consistently by the graph data structure. To enable set :code:`checkMultiEdge` to True. Note that this increases the runtime of the function by O(max(deg(u), deg(v))).

	 	Parameters
	 	----------
		u : int
			Endpoint of edge.
		v : int
			Endpoint of edge.
		w : float, optional
			Edge weight.
		addMissing : bool, optional
			Add missing endpoints if necessary (i.e., increase numberOfNodes). Default: False
		checkMultiEdge : bool, optional
			Check if edge is already present in the graph. If detected, do not insert the edge. Default: False

		Returns
		-------
		bool
			Indicates whether the edge has been added. Is `False` in case :code:`checkMultiEdge` is set to `True` and the new edge would have been a multi-edge.
		"""
		if not (self._this.hasNode(u) and self._this.hasNode(v)):
			if not addMissing:
				raise RuntimeError("Cannot create edge ({0}, {1}) as at least one end point does not exist".format(u,v))

			k = max(u, v)
			if k >= self._this.upperNodeIdBound():
				self._this.addNodes(k - self._this.upperNodeIdBound() + 1)

			if not self._this.hasNode(u):
				self._this.restoreNode(u)

			if not self._this.hasNode(v):
				self._this.restoreNode(v)

		return self._this.addEdge(u, v, w, checkMultiEdge)

	def setWeight(self, u, v, w):
		""" 
		setWeight(u, v, w)
		
		Set the weight of an edge. If the edge does not exist, it will be inserted.

		Parameters
		----------
		u : int
			Endpoint of edge.
		v : int
			Endpoint of edge.
		w : float
			Edge weight.
		"""
		self._this.setWeight(u, v, w)
		return self

	def increaseWeight(self, u, v, w):
		""" 
		increaseWeight(u, v, w)
		
		Increase the weight of an edge. If the edge does not exist, it will be inserted.

		Parameters
		----------
		u : int
			Endpoint of edge.
		v : int
			Endpoint of edge.
		w : float
			Edge weight.
		"""
		self._this.increaseWeight(u, v, w)
		return self

	def removeEdge(self, u, v):
		""" 
		removeEdge(u, v)
		
		Removes the undirected edge {`u`,`v`}.

		Parameters
		----------
		u : int
			Endpoint of edge.
		v : int
			Endpoint of edge.
		"""
		self._this.removeEdge(u, v)
		return self

	def removeAllEdges(self):
		"""
		removeAllEdges()

		Removes all the edges in the graph.
		"""
		self._this.removeAllEdges()

	def removeSelfLoops(self):
		"""
		removeSelfLoops()
		
		Removes all self-loops from the graph.
		"""
		self._this.removeSelfLoops()

	def removeMultiEdges(self):
		""" Removes all multi-edges from the graph.
		"""
		self._this.removeMultiEdges()

	def swapEdge(self, node s1, node t1, node s2, node t2):
		"""
		swapEdge(s1, t1, s2, t2)

		Changes the edge (s1, t1) into (s1, t2) and the edge (s2, t2) into (s2, t1).

		If there are edge weights or edge ids, they are preserved. 
		
		Note
		----
		No check is performed if the swap is actually possible, i.e. does not generate duplicate edges.

		Parameters
		----------
		s1 : int
			Source node of the first edge.
		t1 : int
			Target node of the first edge.
		s2 : int
			Source node of the second edge.
		t2 : int
			Target node of the second edge.
		"""
		self._this.swapEdge(s1, t1, s2, t2)
		return self

	def compactEdges(self):
		"""
		compactEdges()
		
		Compact the edge storage, this should be called after executing many edge deletions.
		"""
		self._this.compactEdges()

	def sortEdges(self):
		"""
		sortEdges()

		Sorts the adjacency arrays by node id. While the running time is linear this
		temporarily duplicates the memory.
		"""
		self._this.sortEdges()

	def hasEdge(self, u, v):
		"""
		hasEdge(u, v) 
		
		Checks if undirected edge {`u`,`v`} exists in the graph.

		Parameters
		----------
		u : int
			Endpoint of edge.
		v : int
			Endpoint of edge.

		Returns
		-------
		bool
			True if the edge exists, False otherwise.
		"""
		return self._this.hasEdge(u, v)

	def weight(self, u, v):
		""" 
		weight(u, v)

		Get edge weight of edge {`u` , `v`}. Returns 0 if edge does not exist.

		Parameters
		----------
		u : int
			Endpoint of edge.
		v : int
			Endpoint of edge.

		Returns
		-------
		float
			Edge weight of edge {`u` , `v`} or 0 if edge does not exist.
		"""
		return self._this.weight(u, v)

	def forNodes(self, object callback):
		""" 
		forNodes(callback)
		
		Experimental node iterator interface

		Parameters
		----------
		callback : object
			Any callable object that takes the parameter node.
		"""
		cdef NodeCallbackWrapper* wrapper
		try:
			wrapper = new NodeCallbackWrapper(callback)
			self._this.forNodes[NodeCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forNodesInRandomOrder(self, object callback):
		""" 
		forNodesInRandomOrder(callback)
		
		Experimental node iterator interface

		Parameters:
		-----------
		callback : object
			Any callable object that takes the parameter node.
		"""
		cdef NodeCallbackWrapper* wrapper
		try:
			wrapper = new NodeCallbackWrapper(callback)
			self._this.forNodesInRandomOrder[NodeCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forNodePairs(self, object callback):
		""" 
		forNodePairs(callback)
		
		Experimental node pair iterator interface

		Parameters
		----------
		callback : object
			Any callable object that takes the parameters tuple(int, int).
			Parameter list refering to (node id, node id).
		"""
		cdef NodePairCallbackWrapper* wrapper
		try:
			wrapper = new NodePairCallbackWrapper(callback)
			self._this.forNodePairs[NodePairCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forEdges(self, object callback):
		""" 
		forEdges(callback)

		Experimental edge iterator interface

		Parameters
		----------
		callback : object
			Any callable object that takes the parameter tuple(int, int, float, int). 
			Parameter list refering to (node id, node id, edge weight, edge id).
		"""
		cdef EdgeCallBackWrapper* wrapper
		try:
			wrapper = new EdgeCallBackWrapper(callback)
			self._this.forEdges[EdgeCallBackWrapper](dereference(wrapper))
		finally:
			del wrapper

	def forEdgesOf(self, node u, object callback):
		""" 
		forEdgesOf(u, callback)
		
		Experimental incident (outgoing) edge iterator interface

		Parameters
		----------
		u : int
			The node of which incident edges shall be passed to the callback
		callback : object
			Any callable object that takes the parameter tuple(int, int, float, int).
			Parameter list refering to (node id, node id, edge weight, edge id).
		"""
		cdef EdgeCallBackWrapper* wrapper
		try:
			wrapper = new EdgeCallBackWrapper(callback)
			self._this.forEdgesOf[EdgeCallBackWrapper](u, dereference(wrapper))
		finally:
			del wrapper

	def forInEdgesOf(self, node u, object callback):
		""" 
		forInEdgesOf(u, callback)
		
		Experimental incident edge iterator interface

		Parameters
		----------
		u : int
			The node of which incident edges shall be passed to the callback
		callback : object
			Any callable object that takes the parameter tuple(int, int, float, int).
			Parameter list refering to (node id, node id, edge weight, edge id).
		"""
		cdef EdgeCallBackWrapper* wrapper
		try:
			wrapper = new EdgeCallBackWrapper(callback)
			self._this.forInEdgesOf[EdgeCallBackWrapper](u, dereference(wrapper))
		finally:
			del wrapper

	def isWeighted(self):
		"""
		isWeighted()

		Returns whether a graph is weighted.

		Returns
		-------
		bool
			True if this graph supports edge weights other than 1.0.
		"""
		return self._this.isWeighted()

	def isDirected(self):
		"""
		isDirected()
		
		Returns whether a graph is directed.
		
		Returns
		-------
		bool
			True if graph is directed.
		"""
		return self._this.isDirected()

	def totalEdgeWeight(self):
		""" 
		totalEdgeWeight()
		
		Get the sum of all edge weights.

		Returns
		-------
		float
			The sum of all edge weights.
		"""
		return self._this.totalEdgeWeight()

	def numberOfSelfLoops(self):
		"""
		numberOfSelfLoops()
		
		Get number of self-loops, i.e. edges {v, v}.

		Returns
		-------
		int
			Number of self-loops.
		"""
		return self._this.numberOfSelfLoops()

	def checkConsistency(self):
		"""
		checkConsistency()

		Check for invalid graph states, such as multi-edges.

		Returns
		-------
		bool
			True if graph contains invalid graph states.
		"""
		return self._this.checkConsistency()

	def iterNodes(self):
		"""
		iterNodes()

		Iterates over the nodes of the graph.
		"""
		it = self._this.nodeRange().begin()
		while it != self._this.nodeRange().end():
			yield dereference(it)
			preincrement(it)

	def iterEdges(self):
		"""
		iterEdges()

		Iterates over the edges of the graph.
		
		For each node u in the graph in ascending node id order,
		the iterator yields the out-edges of u in directed graphs
		and the edges (u,v) in which u < v for undirected graphs.
		
		It does not follow the order of edge ids (if present).
		"""
		it = self._this.edgeRange().begin()
		while it != self._this.edgeRange().end():
			yield dereference(it).u, dereference(it).v
			preincrement(it)

	def iterEdgesWeights(self):
		"""
		iterEdgeWeights()

		Iterates over the edges of the graph and their weights.
		"""
		it = self._this.edgeWeightRange().begin()
		while it != self._this.edgeWeightRange().end():
			yield dereference(it).u, dereference(it).v, dereference(it).weight
			preincrement(it)

	def iterNeighbors(self, u):
		"""
		iterNeighbors(u)

		Iterates over a range of the neighbors of a node.

		Parameters
		----------
		u : int
			The input node.
		"""
		it = self._this.neighborRange(u).begin()
		while it != self._this.neighborRange(u).end():
			yield dereference(it)
			preincrement(it)

	def iterInNeighbors(self, u):
		"""
		iterInNeighbors(u)

		Iterates over a range of the in-neighbors of a node.

		Parameters
		----------
		u : int
			The input node.
		"""
		it = self._this.inNeighborRange(u).begin()
		while it != self._this.inNeighborRange(u).end():
			yield dereference(it)
			preincrement(it)

	def iterNeighborsWeights(self, u):
		"""
		iterNeighborsWeights(u)

		Iterates over a range of the neighbors of a node including the edge weights.
		The iterator is not safe to use with unweighted graphs. To avoid unsafe behavior
		a runtime error will be thrown.

		Parameters
		----------
		u : int
			The input node.
		"""
		if not self._this.isWeighted():
			raise RuntimeError("iterNeighborsWeights: Use this iterator only on weighted graphs.")
		
		it = self._this.weightNeighborRange(u).begin()
		while it != self._this.weightNeighborRange(u).end():
			yield dereference(it)
			preincrement(it)		

	def iterInNeighborsWeights(self, u):
		"""
		iterInNeighborsWeights(u)

		Iterates over a range of the in-neighbors of a node including the edge weights.
		The iterator is not safe to use with unweighted graphs. To avoid unsafe behavior
		a runtime error will be thrown.

		Parameters
		----------
		u : int
			The input node.
		"""
		if not self._this.isWeighted():
			raise RuntimeError("iterInNeighborsWeights: Use this iterator only on weighted graphs.")

		it = self._this.weightInNeighborRange(u).begin()
		while it != self._this.weightInNeighborRange(u).end():
			yield dereference(it)
			preincrement(it)

	def attachNodeAttribute(self, name, ofType):
		"""
		attachNodeAttribute(name, ofType)

		Attaches a node attribute to the graph and returns it.

		.. code-block::
			
			A = G.attachNodeAttribute("attributeIdentifier", ofType)
		
		All values are initially undefined for existing nodes values can be set/get
		by 
		
		.. code-block:: 
		
			A[node] = value # set
			value = A[node] # get

		Getting undefined values raises a ValueError removing a node makes all
		its attributes undefined

		Notes
		-----
		Using node attributes is in experimental state. The API may change in future updates.

		Parameters
		----------
		name   : str
			Name for this attribute
		ofType : type
			Type of the attribute (either int, float, or str)

		Returns
		-------
		networkit.graph.NodeAttribute
			The resulting node attribute container.
		"""
		if not isinstance(name, str):
			raise Exception("Attribute name has to be a string")

		if ofType == int:
			return NodeAttribute(NodeIntAttribute().setThis(self._this.attachNodeIntAttribute(stdstring(name)), &self._this), int)
		elif ofType == float:
			return NodeAttribute(NodeDoubleAttribute().setThis(self._this.attachNodeDoubleAttribute(stdstring(name)), &self._this), float)
		elif ofType == str:
			return NodeAttribute(NodeStringAttribute().setThis(self._this.attachNodeStringAttribute(stdstring(name)), &self._this), str)

	def detachNodeAttribute(self, name):
		"""
		detachNodeAttribute(name)

		Detaches a node attribute from the graph.

		Notes
		-----
		Using node attributes is in experimental state. The API may change in future updates.

		Parameters
		----------
		name : str
			The distinguished name for the attribute to detach.
		"""
		if not isinstance(name, str):
			raise Exception("Attribute name has to be a string")
		self._this.detachNodeAttribute(stdstring(name))

# The following 3 classes NodeIntAttribute, NodeDoubleAttribute and 
# NodeStringAttribute are helper classes which cannot be generalized because
# they map to different C++ classes even if these are generated from the same
# C++ template - this results in some unpleasant code duplication.
# The generic (pure python) wrapper class for the user is NodeAttribute

cdef class NodeIntAttribute:

	cdef setThis(self, _NodeIntAttribute& other, _Graph* G):
		self._this.swap(other)
		self._G = G
		return self

	def __getitem__(self, node):
		try:
			value = self._this.get(node)
		except Exception as e:
			raise ValueError(str(e))
		return value

	def __setitem__(self, node, value):
		try:
			self._this.set(node, value)
		except Exception as e:
			raise ValueError(str(e))

	def __iter__(self):
		try:
			self._iter = self._this.begin()
		except Exception as e:
			raise ValueError(str(e))

		self._stopiter = self._this.end()
		return self

	def __next__(self):
		if self._iter == self._stopiter:
			raise StopIteration()
		val = dereference(self._iter)
		preincrement(self._iter)
		return val


cdef class NodeDoubleAttribute:
	cdef setThis(self, _NodeDoubleAttribute& other, _Graph* G):
		self._this.swap(other)
		self._G = G
		return self

	def __getitem__(self, node):
		try:
			value = self._this.get(node)
		except Exception as e:
			raise ValueError(str(e))
		return value

	def __setitem__(self, node, value):
		try:
			self._this.set(node, value)
		except Exception as e:
			raise ValueError(str(e))

	def __iter__(self):
		try:
			self._iter = self._this.begin()
		except Exception as e:
			raise ValueError(str(e))
		self._stopiter = self._this.end()
		return self

	def __next__(self):
		if self._iter == self._stopiter:
			raise StopIteration()
		val = dereference(self._iter)
		preincrement(self._iter)
		return val

cdef class NodeStringAttribute:

	cdef setThis(self, _NodeStringAttribute& other, _Graph* G):
		self._this.swap(other)
		self._G = G
		return self

	def __getitem__(self, node):
		try:
			value = pystring(self._this.get(node))
		except Exception as e:
			raise ValueError(str(e))
		return value

	def __setitem__(self, node, value):
		try:
			self._this.set(node, stdstring(value))
		except Exception as e:
			raise ValueError(str(e))

	def __iter__(self):
		try:
			self._iter = self._this.begin()
		except Exception as e:
			raise ValueError(str(e))
		self._stopiter = self._this.end()
		return self

	def __next__(self):
		if self._iter == self._stopiter:
			raise StopIteration()
		val = dereference(self._iter)
		val = (val[0], pystring(val[1]))
		preincrement(self._iter)
		return val

class NodeAttribute:
	"""
	Generic class for node attributes returned by networkit.graph.attachNodeAttribute().
	Example of attaching an int attribute to a graph g:

	.. code-block::

		att = g.attachNodeAttribute("name", int)`

	Set/get attributes of a single node 'u' with the [] operator:

	.. code-block::

		att[u] = 0
		att_val = att[u] # 'att_val' is 0

	Iterate over all the values of an attribute:

	.. code-block::

		for u, val in att:
			# The attribute value of node `u` is `val`.

	Notes
	-----
	Using node attributes is in experimental state. The API may change in future updates.
	"""

	def __init__(self, typedNodeAttribute, type):
		self.attr = typedNodeAttribute
		self.type = type

	def __getitem__(self, node):
		return self.attr[node]

	def __setitem__(self, index, value):
		if not isinstance(value, self.type):
			raise Exception("Wrong Attribute type")
		self.attr[index] = value

	def __iter__(self):
		self._iter = iter(self.attr)
		return self

	def __next__(self):
		return next(self._iter)

cdef cppclass EdgeCallBackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(node u, node v, edgeweight w, edgeid eid):
		cdef bool_t error = False
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
		cdef bool_t error = False
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
		cdef bool_t error = False
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
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(u, v)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef class SpanningForest:
	""" 
	SpanningForest(G, nodes)
	
	Generates a spanning forest for a given graph

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	nodes : list(int)
		A subset of nodes of `G` which induce the subgraph.
	"""
	cdef _SpanningForest* _this
	cdef Graph _G

	def __cinit__(self, Graph G not None):
		self._G = G
		self._this = new _SpanningForest(G._this)


	def __dealloc__(self):
		del self._this

	def run(self):
		"""
		run()

		Executes the algorithm.
		"""
		self._this.run()
		return self

	def getForest(self):
		"""
		getForest()

		Returns the spanning forest.

		Returns
		-------
		networkit.Graph
			The computed spanning forest.
		"""
		return Graph().setThis(self._this.getForest())

cdef class RandomMaximumSpanningForest(Algorithm):
	"""
	RandomMaximumSpanningForest(G, attributes)

	Computes a random maximum-weight spanning forest using Kruskal's algorithm by randomizing the order of edges of the same weight.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	attribute : list(int) or list(float)
		If given, this edge attribute is used instead of the edge weights.
	"""

	def __cinit__(self, Graph G not None, vector[double] attribute = vector[double]()):
		self._G = G
		if attribute.empty():
			self._this = new _RandomMaximumSpanningForest(G._this)
		else:
			self._attribute = move(attribute)
			self._this = new _RandomMaximumSpanningForest(G._this, self._attribute)

	def getMSF(self, bool_t move):
		"""
		getMSF(move)

		Gets the calculated maximum-weight spanning forest as graph.

		Parameters
		----------
		move : bool
			If the graph shall be moved out of the algorithm instance.

		Returns
		-------
		networkit.Graph
			The calculated maximum-weight spanning forest.
		"""
		return Graph().setThis((<_RandomMaximumSpanningForest*>(self._this)).getMSF(move))

	def getAttribute(self, bool_t move = False):
		"""
		getAttribute(move=False)

		Get a bool attribute that indicates for each edge if it is part of the calculated maximum-weight spanning forest.
		This attribute is only calculated and can thus only be request if the supplied graph has edge ids.

		Parameters
		----------
		move : bool, optional
			If the attribute shall be moved out of the algorithm instance. Default: False

		Returns
		-------
		list(bool)
			The list with the bool attribute for each edge.
		"""
		return (<_RandomMaximumSpanningForest*>(self._this)).getAttribute(move)

	def inMSF(self, node u, node v = _none):
		"""
		inMSF(u, v = None)

		Checks if the edge (u, v) or the edge with id u is part of the calculated maximum-weight spanning forest.

		Parameters
		----------
		u : int
			The first node of the edge to check or the edge id of the edge to check.
		v : int, optional
			The second node of the edge to check (only if u is not an edge id). Default: None

		Returns
		-------
		bool
			If the edge is part of the calculated maximum-weight spanning forest.
		"""
		if v == _none:
			return (<_RandomMaximumSpanningForest*>(self._this)).inMSF(u)
		else:
			return (<_RandomMaximumSpanningForest*>(self._this)).inMSF(u, v)

cdef class UnionMaximumSpanningForest(Algorithm):
	"""
	UnionMaximumSpanningForest(G, attribute)

	Union maximum-weight spanning forest algorithm, computes the union of all maximum-weight spanning forests using Kruskal's algorithm.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	attribute : list(int) or list(float)
		If given, this edge attribute is used instead of the edge weights.
	"""

	def __cinit__(self, Graph G not None, vector[double] attribute = vector[double]()):
		self._G = G

		if attribute.empty():
			self._this = new _UnionMaximumSpanningForest(G._this)
		else:
			self._this = new _UnionMaximumSpanningForest(G._this, attribute)

	def getUMSF(self, bool_t move = False):
		"""
		getUMSF(move=False)

		Gets the union of all maximum-weight spanning forests as graph.

		Parameters
		----------
		move : bool, optional
			If the graph shall be moved out of the algorithm instance. Default: False

		Returns
		-------
		networkit.Graph
			The calculated union of all maximum-weight spanning forests.
		"""
		return Graph().setThis((<_UnionMaximumSpanningForest*>(self._this)).getUMSF(move))

	def getAttribute(self, bool_t move = False):
		"""
		getAttribute(move=False)
		
		Get a bool attribute that indicates for each edge if it is part of any maximum-weight spanning forest.

		This attribute is only calculated and can thus only be request if the supplied graph has edge ids.

		Parameters
		----------
		move : bool, optional
			If the attribute shall be moved out of the algorithm instance. Default: False

		Returns
		-------
		list(bool)
			The list with the bool attribute for each edge.
		"""
		return (<_UnionMaximumSpanningForest*>(self._this)).getAttribute(move)

	def inUMST(self, node u, node v = _none):
		"""
		inUMST(u, v=None)

		Checks if the edge (u, v) or the edge with id u is part of any maximum-weight spanning forest.

		Parameters
		----------
		u : int
			The first node of the edge to check or the edge id of the edge to check.
		v : int, optional
			The second node of the edge to check (only if u is not an edge id). Default: None

		Returns
		-------
		bool
			If the edge is part of any maximum-weight spanning forest.
		"""
		if v == _none:
			return (<_UnionMaximumSpanningForest*>(self._this)).inUMSF(u)
		else:
			return (<_UnionMaximumSpanningForest*>(self._this)).inUMSF(u, v)
