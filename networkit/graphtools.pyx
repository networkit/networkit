# distutils: language=c++

from libc.stdint cimport uint64_t
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node
ctypedef double edgeweight

from .graph cimport _Graph, Graph

cdef extern from "<networkit/graph/GraphTools.hpp>" namespace "NetworKit::GraphTools":

	count maxDegree(_Graph G) nogil except +
	count maxInDegree(_Graph G) nogil except +
	edgeweight maxWeightedDegree(_Graph G) nogil except +
	edgeweight maxWeightedInDegree(_Graph G) nogil except +
	node randomNode(_Graph G) nogil except +
	node randomNeighbor(_Graph G, node u) nogil except +
	pair[node, node] randomEdge(_Graph G, bool_t uniformDistribution) nogil except +
	vector[pair[node, node]] randomEdges(_Graph G, count numEdges) nogil except +
	pair[count, count] size(_Graph G) nogil except +
	double density(_Graph G) nogil except +
	double volume(_Graph G) nogil except +
	double volume[InputIt](_Graph G, InputIt first, InputIt last) nogil except +
	double inVolume[InputIt](_Graph G, InputIt first, InputIt last) nogil except +
	_Graph copyNodes(_Graph G) nogil except +
	_Graph toUndirected(_Graph G) nogil except +
	_Graph toUnweighted(_Graph G) nogil except +
	_Graph toWeighted(_Graph G) nogil except +
	_Graph subgraphFromNodes(_Graph G, unordered_set[node], bool_t, bool_t) nogil except +
	void append(_Graph G, _Graph G1) nogil except +
	void merge(_Graph G, _Graph G1) nogil except +
	void removeEdgesFromIsolatedSet[InputIt](_Graph G, InputIt first, InputIt last) except +
	_Graph getCompactedGraph(_Graph G, unordered_map[node,node]) nogil except +
	_Graph transpose(_Graph G) nogil except +
	unordered_map[node,node] getContinuousNodeIds(_Graph G) nogil except +
	unordered_map[node,node] getRandomContinuousNodeIds(_Graph G) nogil except +

cdef class GraphTools:

	@staticmethod
	def maxDegree(Graph G):
		"""
		Returns the maximum out-degree of the graph.

		Parameters:
		-----------
		G : networkit.Graph
			The input graph.

		Returns:
		--------
		count
			The maximum out-degree of the graph.
		"""
		return maxDegree(G._this)

	@staticmethod
	def maxInDegree(Graph G):
		"""
		Returns the maximum in-degree of the graph.

		Parameters:
		-----------
		G : networkit.Graph
			The input graph.

		Returns:
		--------
		count
			The maximum in-degree of the graph.
		"""
		return maxInDegree(G._this)

	@staticmethod
	def maxWeightedDegree(Graph G):
		"""
		Returns the maximum weighted out-degree of the graph.

		Parameters:
		-----------
		G : networkit.Graph
			The input graph.

		Returns:
		--------
		edgeweight
			The maximum weighted out-degree of the graph.
		"""
		return maxWeightedDegree(G._this)

	@staticmethod
	def maxWeightedInDegree(Graph G):
		"""
		Returns the maximum weighted in-degree of the graph.

		Parameters:
		-----------
		G : networkit.Graph
			The input graph.

		Returns:
		--------
		edgeweight
			The maximum weighted in-degree of the graph.
		"""
		return maxWeightedInDegree(G._this)

	@staticmethod
	def randomNode(Graph G):
		"""
		Returns a random node of the input graph.

		Parameters:
		-----------
		G : networkit.Graph
			The input graph.

		Returns:
		--------
		node
			A random node.
		"""
		return randomNode(G._this)

	@staticmethod
	def randomNeighbor(Graph G, node u):
		"""
		Returns a random neighbor of node `u`.

		Parameters:
		-----------
		G : networkit.Graph
			The input graph.
		u : node
			A node in `G`.

		Returns:
		--------
		node
			A random neighbor of `u`.
		"""
		return randomNeighbor(G._this, u)

	@staticmethod
	def randomEdge(Graph G, uniformDistribution = False):
		""" Get a random edge of the graph.

		Parameters:
		-----------
		G : networkit.Graph
			The input graph.
		uniformDistribution : bool
			If the distribution of the edge shall be uniform.

		Returns:
		--------
		pair
			Random edge.

		Notes
		-----
		Fast, but not uniformly random if uniformDistribution is not set,
		slow and uniformly random otherwise.
		"""
		return randomEdge(G._this, uniformDistribution)

	@staticmethod
	def randomEdges(Graph G, numEdges):
		"""
		Returns a list with numEdges random edges. The edges are chosen uniformly at random.

		Parameters:
		-----------
		G : networkit.Graph
			The input graph.
		numEdges : count
			The number of edges to choose.

		Returns:
		--------
		list of pairs
			List of with `numEdges` random edges.
		"""
		return randomEdges(G._this, numEdges)

	@staticmethod
	def append(Graph G, Graph G1):
		"""
		Appends graph `G1` to graph `G` as a new subgraph. Performs node id remapping.

		Parameters:
		-----------
		G : networkit.Graph
			Graph where `G1` will be appended to.
		G1 : networkit.Graph
			Graph that will be appended to `G`.
		"""
		append(G._this, G1._this)

	@staticmethod
	def merge(Graph G, Graph G1):
		"""
		Modifies graph `G` to be the union of it and graph `G1`.
		Nodes with the same ids are identified with each other.

		Parameters:
		-----------
		G : networkit.Graph
			Result of the merge.
		G1 : networkit.Graph
			Graph that will be merged with `G`.
		"""
		merge(G._this, G1._this)

	@staticmethod
	def removeEdgesFromIsolatedSet(Graph graph, nodes):
		"""
		Efficiently removes all the edges adjacent to a set of nodes that is
		not connected to the rest of the graph. This is meant to optimize the
		Kadabra algorithm.

		Parameters:
		-----------
		G : networkit.Graph
			The input graph.
		nodes : list
			Isolates set of nodes from where the edges will be removed.
		"""
		cdef vector[node] isolatedSet

		try:
			isolatedSet = <vector[node]?>nodes
		except TypeError:
			raise RuntimeError("Error, nodes must be a list of nodes.")
		removeEdgesFromIsolatedSet[vector[node].iterator](graph._this,
				isolatedSet.begin(), isolatedSet.end())

	@staticmethod
	def toUndirected(Graph graph):
		"""
		Returns an undirected copy of the input graph.

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.

		Returns:
		--------
		graph : networkit.Graph
			Undirected copy of the input graph.
		"""
		return Graph().setThis(toUndirected(graph._this))

	@staticmethod
	def toUnweighted(Graph graph):
		"""
		Returns an unweighted copy of the input graph.

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.

		Returns:
		--------
		graph : networkit.Graph
			Unweighted copy of the input graph.
		"""
		return Graph().setThis(toUnweighted(graph._this))

	@staticmethod
	def toWeighted(Graph graph):
		"""
		Returns a weighted copy of the input graph.

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.

		Returns:
		--------
		graph : networkit.Graph
			Weighted copy of the input graph.
		"""
		return Graph().setThis(toWeighted(graph._this))

	@staticmethod
	def size(Graph graph):
		"""
		Return the size of the graph.

		Returns:
		--------
		tuple
			a pair (n, m) where n is the number of nodes and m is the number of edges.
		"""
		return size(graph._this)

	@staticmethod
	def density(Graph graph):
		"""
		Get the density of the input graph.

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.

		Returns:
		--------
		double
			The density of the input graph.
		"""
		return density(graph._this)

	@staticmethod
	def volume(Graph graph, nodes = None):
		"""
		Get the volume (for all outgoing edges) of a graph. If a list of nodes of the graph
		is given, the volume for the corresponding subgraph is computed.

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.
		nodes : list
			(Optional) List of nodes from the graph.

		Returns:
		--------
		double
			The volume of the subgraph.
		"""

		cdef vector[node] cNodes

		if nodes is not None:
			try:
				cNodes = <vector[node]?>nodes
				return volume[vector[node].iterator](graph._this, cNodes.begin(), cNodes.end())
			except TypeError:
				raise RuntimeError("Error, nodes must be a list of nodes.")
		else:
			return volume(graph._this)

	@staticmethod
	def inVolume(Graph graph, nodes):
		"""
		Get the inVolume (for all incoming edges) of a subgraph, defined by the 
		input graph and a corresponding subset of nodes.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		nodes : vector[node]
			A vector of nodes from the graph. 

		Returns
		-------
		double
			The inVolume of the input graph.
		"""

		cdef vector[node] cNodes

		try:
			cNodes = <vector[node]?>nodes
			return inVolume[vector[node].iterator](graph._this, cNodes.begin(), cNodes.end())
		except TypeError:
			raise RuntimeError("Error, nodes must be a list of nodes.")

	@staticmethod
	def copyNodes(Graph graph):
		"""
		Copies all nodes of the input graph to a new graph (edges are not copied).

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.

		Returns:
		--------
		graph : networkit.Graph
			Graph with the same nodes as the input graph (and without any edge).
		"""
		return Graph().setThis(copyNodes(graph._this))

	@staticmethod
	def subgraphFromNodes(Graph graph, nodes, includeOutNeighbors=False, includeInNeighbors=False):
		"""
		Returns an induced subgraph of the input graph (including potential edge
		weights/directions).

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.
		nodes : set
			Nodes in the induced subgraph.
		includeOutNeighbors : bool
			If set to true, out-neighbors will also be included.
		includeInNeighbors : bool
			If set to true, in-neighbors will also be included.

		Returns:
		--------
		graph : networkit.Graph
			Induced subgraph.
		"""
		return Graph().setThis(subgraphFromNodes(
			graph._this, nodes, includeOutNeighbors, includeInNeighbors))

	@staticmethod
	def transpose(Graph graph):
		"""
		Returns the transpose of the input graph. The graph must be directed.

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.

		Returns:
		graph : networkit.Graph
			Transpose of the input graph.
		"""
		return Graph().setThis(transpose(graph._this))

	@staticmethod
	def getCompactedGraph(Graph graph, nodeIdMap):
		"""
		Computes a graph with the same structure but with continuous node ids.

		Parameters:
		-----------
		graph : networkit.Graph
			The graph to be compacted.
		nodeIdMap:
			The map providing the information about the node ids.

		Returns:
		--------
		networkit.Graph
			The compacted graph
		"""
		cdef unordered_map[node,node] cNodeIdMap
		for key in nodeIdMap:
			cNodeIdMap[key] = nodeIdMap[key]
		return Graph().setThis(getCompactedGraph(graph._this,cNodeIdMap))

	@staticmethod
	def getContinuousNodeIds(Graph graph):
		"""
		Computes a map of node ids to continuous node ids.

		Parameters:
		-----------
		graph : networkit.Graph
			The graph of which the node id map is wanted.
		Returns:
		--------
			Returns the node id map
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

		Parameters:
		-----------
		graph : networkit.Graph
			The graph of which the node id map is wanted.
		Returns:
		--------
			Returns the node id map
		"""
		cdef unordered_map[node,node] cResult
		with nogil:
			cResult = getRandomContinuousNodeIds(graph._this)
		result = dict()
		for elem in cResult:
			result[elem.first] = elem.second
		return result
