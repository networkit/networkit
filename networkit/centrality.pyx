cdef extern from "<networkit/centrality/DegreeCentrality.hpp>":

	cdef cppclass _DegreeCentrality "NetworKit::DegreeCentrality" (_Centrality):
		_DegreeCentrality(_Graph, bool_t normalized, bool_t outdeg, bool_t ignoreSelfLoops) except +

cdef class Centrality(Algorithm):
	""" Abstract base class for centrality measures"""


	def __init__(self, *args, **kwargs):
		if type(self) == Centrality:
			raise RuntimeError("Error, you may not use Centrality directly, use a sub-class instead")

	def __dealloc__(self):
		self._G = None # just to be sure the graph is deleted

	def scores(self):
		"""
		Returns
		-------
		list
			the list of all scores
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).scores()

	def score(self, v):
		"""
		Returns
		-------
		the score of node v
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).score(v)

	def ranking(self):
		"""
		Returns
		-------
		dictionary
			a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).ranking()

	def maximum(self):
		"""
		Returns
		-------
		the maximum theoretical centrality score for the given graph
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).maximum()

	def centralization(self):
		"""
		Compute the centralization of a network with respect to some centrality measure.
		The centralization of any network is a measure of how central its most central
		node is in relation to how central all the other nodes are.
		Centralization measures then (a) calculate the sum in differences
		in centrality between the most central node in a network and all other nodes;
		and (b) divide this quantity by the theoretically largest such sum of
		differences in any network of the same size.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).centralization()

cdef class DegreeCentrality(Centrality):
	""" Node centrality index which ranks nodes by their degree.
	Optional normalization by maximum degree. The run() method runs in O(m) time, where m is the number of
	edges in the graph.

	DegreeCentrality(G, normalized=False, outDeg=True, ignoreSelfLoops=True)

	Constructs the DegreeCentrality class for the given Graph `G`. If the scores should be normalized,
	then set `normalized` to True.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	normalized : bool, optional
		Normalize centrality values in the interval [0,1].
		outdeg : bool, optional
		If set to true, computes the centrality based on out-degrees, otherwise based on the in-degrees.
		ignoreSelfLoops : bool, optional
		If set to true, self loops will not be taken into account.
	"""

	def __cinit__(self, Graph G, bool_t normalized=False, bool_t outDeg = True, bool_t ignoreSelfLoops=True):
		self._G = G
		self._this = new _DegreeCentrality(G._this, normalized, outDeg, ignoreSelfLoops)

cdef extern from "<networkit/centrality/LocalPartitionCoverage.hpp>":

	cdef cppclass _LocalPartitionCoverage "NetworKit::LocalPartitionCoverage" (_Centrality):
		_LocalPartitionCoverage(_Graph, _Partition) except +

cdef class LocalPartitionCoverage(Centrality):
	"""
	The local partition coverage is the amount of neighbors of a node u that are in the same partition as u.
	The running time of the run() method is O(m), where m is the number of edges in the graph.

	LocalPartitionCoverage(G, P)

	Parameters
	----------
	G : networkit.Graph
		The graph.
	P : networkit.Partition
		The partition to use
	"""
	cdef Partition _P

	def __cinit__(self, Graph G not None, Partition P not None):
		self._G = G
		self._P = P
		self._this = new _LocalPartitionCoverage(G._this, P._this)

