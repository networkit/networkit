# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set

from .graph cimport _Graph, Graph
from .structures cimport count, index, node, edgeweight

cdef extern from "<networkit/graph/GraphTools.hpp>" namespace "NetworKit::GraphTools":

	count maxDegree(_Graph G) nogil except +
	count maxInDegree(_Graph G) nogil except +
	edgeweight maxWeightedDegree(_Graph G) nogil except +
	edgeweight maxWeightedInDegree(_Graph G) nogil except +
	node randomNode(_Graph G) nogil except +
	vector[node] randomNodes(_Graph G, count n) nogil except +
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
	_Graph subgraphFromNodes[InputIt](_Graph G, InputIt first, InputIt last, bool_t compact) nogil except +
	_Graph subgraphAndNeighborsFromNodes(_Graph G, unordered_set[node], bool_t, bool_t) nogil except +
	void append(_Graph G, _Graph G1) nogil except +
	void merge(_Graph G, _Graph G1) nogil except +
	void removeEdgesFromIsolatedSet[InputIt](_Graph G, InputIt first, InputIt last) except +
	_Graph getCompactedGraph(_Graph G, unordered_map[node,node]) nogil except +
	_Graph transpose(_Graph G) nogil except +
	unordered_map[node,node] getContinuousNodeIds(_Graph G) nogil except +
	unordered_map[node,node] getRandomContinuousNodeIds(_Graph G) nogil except +
	void sortEdgesByWeight(_Graph G, bool_t) nogil except +
	vector[node] topologicalSort(_Graph G) nogil except +
	node augmentGraph(_Graph G) nogil except +
	pair[_Graph, node] createAugmentedGraph(_Graph G) nogil except +

cdef class GraphTools:

	@staticmethod
	def maxDegree(Graph G):
		"""
		maxDegree(G)

		Returns the maximum out-degree of the graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.

		Returns
		-------
		int
			The maximum out-degree of the graph.
		"""
		return maxDegree(G._this)

	@staticmethod
	def maxInDegree(Graph G):
		"""
		maxInDegree(G)

		Returns the maximum in-degree of the graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.

		Returns
		-------
		int
			The maximum in-degree of the graph.
		"""
		return maxInDegree(G._this)

	@staticmethod
	def maxWeightedDegree(Graph G):
		"""
		maxWeightedDegree(G)

		Returns the maximum weighted out-degree of the graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.

		Returns
		-------
		float
			The maximum weighted out-degree of the graph.
		"""
		return maxWeightedDegree(G._this)

	@staticmethod
	def maxWeightedInDegree(Graph G):
		"""
		maxWeightedInDegree(G)

		Returns the maximum weighted in-degree of the graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.

		Returns
		-------
		float
			The maximum weighted in-degree of the graph.
		"""
		return maxWeightedInDegree(G._this)

	@staticmethod
	def randomNode(Graph G):
		"""
		randomNode(G)

		Returns a random node of the input graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.

		Returns
		-------
		int
			A random node.
		"""
		return randomNode(G._this)

	@staticmethod
	def randomNodes(Graph G, count n):
		"""
		randomNodes(G, n)

		Returns n distinct random nodes of the input graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		n : int
			The number of desired nodes.

		Returns
		-------
		list(int)
			A list of distinct random nodes.
		"""
		return randomNodes(G._this, n)

	@staticmethod
	def randomNeighbor(Graph G, node u):
		"""
		randomNeighbor(G, u)

		Returns a random neighbor of node `u`.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		u : int
			A node in `G`.

		Returns
		-------
		int
			A random neighbor of `u`.
		"""
		return randomNeighbor(G._this, u)

	@staticmethod
	def randomEdge(Graph G, uniformDistribution = False):
		""" 
		randomEdge(G, uniformDistribution=False)

		Get a random edge of the graph.

		Notes
		-----
		Fast, but not uniformly random if uniformDistribution is not set,
		slow and uniformly random otherwise.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		uniformDistribution : bool, optional
			If the distribution of the edge shall be uniform. Default: False

		Returns
		-------
		tuple(int, int)
			Random edge.
		"""
		return randomEdge(G._this, uniformDistribution)

	@staticmethod
	def randomEdges(Graph G, numEdges):
		"""
		randomEdges(G, numEdges)

		Returns a list with numEdges random edges. The edges are chosen uniformly at random.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		numEdges : int
			The number of edges to choose.

		Returns
		-------
		list(tuple(int, int))
			List of with `numEdges` random edges.
		"""
		return randomEdges(G._this, numEdges)

	@staticmethod
	def append(Graph G, Graph G1):
		"""
		append(G, G1)

		Appends graph `G1` to graph `G` as a new subgraph. Performs node id remapping.

		Parameters
		----------
		G : networkit.Graph
			Graph where `G1` will be appended to.
		G1 : networkit.Graph
			Graph that will be appended to `G`.
		"""
		append(G._this, G1._this)

	@staticmethod
	def merge(Graph G, Graph G1):
		"""
		merge(G, G1)

		Modifies graph `G` to be the union of it and graph `G1`.
		Nodes with the same ids are identified with each other.

		Parameters
		----------
		G : networkit.Graph
			Result of the merge.
		G1 : networkit.Graph
			Graph that will be merged with `G`.
		"""
		merge(G._this, G1._this)

	@staticmethod
	def removeEdgesFromIsolatedSet(Graph graph, nodes):
		"""
		removeEdgesFromIsolatedSet(graph, nodes)

		Efficiently removes all the edges adjacent to a set of nodes that is
		not connected to the rest of the graph. This is meant to optimize the
		Kadabra algorithm.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		nodes : list(int)
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
		toUndirected(graph)

		Returns an undirected copy of the input graph.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.

		Returns
		-------
		graph : networkit.Graph
			Undirected copy of the input graph.
		"""
		return Graph().setThis(toUndirected(graph._this))

	@staticmethod
	def toUnweighted(Graph graph):
		"""
		toUnweighted(graph)

		Returns an unweighted copy of the input graph.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.

		Returns
		-------
		graph : networkit.Graph
			Unweighted copy of the input graph.
		"""
		return Graph().setThis(toUnweighted(graph._this))

	@staticmethod
	def toWeighted(Graph graph):
		"""
		toWeighted(graph)

		Returns a weighted copy of the input graph.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.

		Returns
		-------
		graph : networkit.Graph
			Weighted copy of the input graph.
		"""
		return Graph().setThis(toWeighted(graph._this))

	@staticmethod
	def size(Graph graph):
		"""
		size(graph)
		
		Return the size of the graph.

		Returns
		-------
		tuple(int, int)
			a pair (n, m) where n is the number of nodes and m is the number of edges.
		"""
		return size(graph._this)

	@staticmethod
	def density(Graph graph):
		"""
		density(graph)

		Get the density of the input graph.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.

		Returns
		-------
		float
			The density of the input graph.
		"""
		return density(graph._this)

	@staticmethod
	def volume(Graph graph, nodes = None):
		"""
		volume(graph, nodes = None)
		
		Get the volume (for all outgoing edges) of a graph. If a list of nodes of the graph
		is given, the volume for the corresponding subgraph is computed.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		nodes : list(int), optional
			List of nodes from the graph.

		Returns
		-------
		float
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
		inVolume(graph, nodes)

		Get the inVolume (for all incoming edges) of a subgraph, defined by the
		input graph and a corresponding subset of nodes.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		nodes : list(int)
			A vector of nodes from the graph.

		Returns
		-------
		float
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
		copyNodes(graph)

		Copies all nodes of the input graph to a new graph (edges are not copied).

		Parameters
		----------
		graph : networkit.Graph
			The input graph.

		Returns
		-------
		graph : networkit.Graph
			Graph with the same nodes as the input graph (and without any edge).
		"""
		return Graph().setThis(copyNodes(graph._this))

	@staticmethod
	def subgraphFromNodes(Graph graph, vector[node] nodes, includeOutNeighbors=False, includeInNeighbors=False, bool_t compact = False):
		"""
		subgraphFromNodes(graph, list(int) nodes, includeOutNeighbors=False, includeInNeighbors=False, compact = False)
		
		Returns an induced subgraph of this graph (including potential edge
		weights/directions)

		The subgraph contains all nodes in Nodes  and all edges which
		have one end point in Nodes and the other in Nodes.

		Parameters:
		-----------
		graph : networkit.Graph
			The input graph.
		nodes : list(int)
			Nodes in the induced subgraph.
		includeOutNeighbors : bool, optional
			DEPRECATED. Use subgraphAndNeighborsFromNodes instead.
			If set to true, out-neighbors will also be included. Default: False
		includeInNeighbors : bool, optional
			DEPRECATED. Use subgraphAndNeighborsFromNodes instead.
			If set to true, in-neighbors will also be included. Default: False
		compact : bool, optional
			Indicates whether the resulting graph shall have compact, continuous node ids.
			If False node ids of the input graph are kept. Default: False

		Returns
		-------
		graph : networkit.Graph
			Induced subgraph.
		"""
		# Deprecated compatibility wrapper. We use "vector" to
		# preserve the sorting of the nodes for compact
		# subgraphs and only convert to unordered_set when
		# needed.
		cdef unordered_set[node] nodeSet

		if includeInNeighbors or includeOutNeighbors:
			if compact:
				raise RuntimeError("Compaction is not supported with includeOutNeighbors or includeInNeighbors")
			nodeSet.insert(nodes.begin(), nodes.end())
			return Graph().setThis(subgraphAndNeighborsFromNodes(
			    	graph._this, nodeSet, includeOutNeighbors, includeInNeighbors))
		else:
			return Graph().setThis(subgraphFromNodes(
			    	graph._this, nodes.begin(), nodes.end(), compact))

	@staticmethod
	def subgraphAndNeighborsFromNodes(Graph graph, nodes, includeOutNeighbors=False, includeInNeighbors=False):
		"""
		subgraphAndNeighborsFromNodes(graph, nodes, includeOutNeighbors=False, includeInNeighbors=False)

		Returns an induced subgraph of this graph (including potential edge
		weights/directions)

		There a two relevant sets of nodes:

		- Nodes are such passed as arguments.
		- Neighbors are empty by default.

		The subgraph contains all nodes in Nodes + Neighbors and all edges which
		have one end point in Nodes and the other in Nodes or Neighbors.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		nodes : list(int)
			Nodes in the induced subgraph.
		includeOutNeighbors : bool, optional
			If set to True, out-neighbors will also be included. Default: False
		includeInNeighbors : bool, optional
			If set to True, in-neighbors will also be included. Default: False

		Returns
		-------
		graph : networkit.Graph
			Induced subgraph.
		"""
		return Graph().setThis(subgraphAndNeighborsFromNodes(
			graph._this, nodes, includeOutNeighbors, includeInNeighbors))

	@staticmethod
	def transpose(Graph graph):
		"""
		transpose(graph)

		Returns the transpose of the input graph. The graph must be directed.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.

		Returns
		-------
		graph : networkit.Graph
			Transpose of the input graph.
		"""
		return Graph().setThis(transpose(graph._this))

	@staticmethod
	def getCompactedGraph(Graph graph, nodeIdMap):
		"""
		getCompactedGraph(graph, nodeIdMap)

		Computes a graph with the same structure but with continuous node ids.

		Parameters
		----------
		graph : networkit.Graph
			The graph to be compacted.
		nodeIdMap : list(int)
			The map providing the information about the node ids.

		Returns
		-------
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
		getContinuousNodeIds(graph)

		Computes a map of node ids to continuous node ids.

		Parameters
		----------
		graph : networkit.Graph
			The graph of which the node id map is wanted.

		Returns
		-------
		list(int)
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
		getRandomContinuousNodeIds(graph):

		Computes a map of node ids to continuous, randomly permutated node ids.

		Parameters
		----------
		graph : networkit.Graph
			The graph of which the node id map is wanted.

		Returns
		-------
		list(int)
			Returns the node id map
		"""
		cdef unordered_map[node,node] cResult
		with nogil:
			cResult = getRandomContinuousNodeIds(graph._this)
		result = dict()
		for elem in cResult:
			result[elem.first] = elem.second
		return result

	@staticmethod
	def sortEdgesByWeight(Graph G, decreasing = False):
		"""
		sortEdgesByWeight(G, decreasing = False)

		Sorts the adjacency arrays by edge weight.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		decreasing : bool, optional
			If True adjacency arrays are sorted by non-increasing edge weights, if False
			adjacency arrays are sorted by non-decreasing edge weights. Ties are broken
			by using node ids. Default: False
		"""
		sortEdgesByWeight(G._this, decreasing)

	@staticmethod
	def topologicalSort(Graph G):
		"""
		topologicalSort(G)

		Given a directed graph G, the topology sort algorithm creates one valid topology order of nodes.
		Undirected graphs are not accepted as input, since a topology sort is a linear ordering of vertices 
		such that for every edge u -> v, node u comes before v in the ordering.

		Parameters
		----------
		G : networkit.Graph
			The directed input graph.
		"""
		return topologicalSort(G._this)

	@staticmethod
	def augmentGraph(Graph G):
		"""
		Augments the input graph in-place as required by ForestCentrality. With respect to the input
		graph G, the augmented graph has a new root node connected to all the other nodes in the graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph (undirected).

		Returns
		-------
		int
			Returns the node id of the new root node.
		"""
		return augmentGraph(G._this)

	@staticmethod
	def createAugmentedGraph(Graph G):
		"""
		Constructs an augmented graph as required by ForestCentrality. With respect to the input
		graph G, the augmented graph has a new root node connected to all the other nodes in the
		graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph (undirected).

		Returns
		-------
		tuple(networkit.Graph, int)
			Returns a tuple (G, root) where G is the augmented graph and root is the id of the root
			node.
		"""
		result = createAugmentedGraph(G._this)
		return Graph().setThis(result.first), result.second
