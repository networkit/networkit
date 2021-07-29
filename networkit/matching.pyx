# distutils: language=c++

cdef class Matching:
	""" Implements a graph matching.

		Matching(z=0)

		Create a new matching data structure for `z` elements.

		Parameters:
		-----------
		z : index, optional
			Maximum number of nodes.
	"""
	def __cinit__(self, index z=0):
		self._this = move(_Matching(z))

	cdef setThis(self,  _Matching& other):
		swap[_Matching](self._this,  other)
		return self

	def match(self, node u, node v):
		"""
		Set two nodes u and v as each others matching partners.

		Parameters:
		-----------
		u : node
		v : node
		"""
		self._this.match(u,v)

	def unmatch(self, node u,  node v):
		"""
		Reset the two nodes u and v to unmatched.

		Parameters:
		-----------
		u : node
		v : node
		"""
		self._this.unmatch(u, v)

	def isMatched(self, node u):
		"""
		Check if node u is matched.

		Parameters:
		-----------
		u : node

		Returns:
		-----------
		bool
			true if u has a matching partner.
		"""
		return self._this.isMatched(u)

	def areMatched(self, node u, node v):
		"""
		Check if the nodes u and v are matched to each other.

		Parameters:
		-----------
		u : node
		v : node

		Returns:
		-----------
		bool
			true if node u and v are matched to each other.
		"""
		return self._this.areMatched(u,v)

	def isProper(self, Graph G):
		"""
		Check whether this is a proper matching in the graph, i.e. there is an edge between
		each pair and if the matching partner of v is u then the matching partner of u is v.

		Parameters:
		-----------
		G : networkit.Graph
			Graph to be checked.
		Returns:
		-----------
		bool
			true if this is a proper matching.
		"""
		return self._this.isProper(G._this)

	def size(self, Graph G):
		"""
		Get the number of edges in this matching.

		Parameters:
		-----------
		G : networkit.Graph

		Returns:
		-----------
		int
			total number of edges in the matching.
		"""
		return self._this.size(G._this)

	def mate(self, node v):
		"""
		Get the matched neighbor of v if it exists, otherwise +inf.

		Parameters:
		-----------
		v : node

		Returns:
		-----------
		int
			matching partner of v if it exists, otherwise +inf.
		"""
		return self._this.mate(v)

	def weight(self, Graph G):
		"""
		Get total weight of edges in this matching.

		Parameters:
		-----------
		G : networkit.Graph
			The corresponding graph.

		Returns:
		-----------
		float
			Total weight of edges in this matching.
		"""
		return self._this.weight(G._this)

	def toPartition(self, Graph G):
		"""
		Convert matching to a Partition object where matched nodes belong to the same subset
		and unmatched nodes belong to a singleton subset.

		Parameters:
		-----------
		G : networkit.Graph
			The corresponding graph.

		Returns:
		-----------
		networkit.Partition

		"""
		return Partition().setThis(self._this.toPartition(G._this))

	def getVector(self):
		""" Get the vector storing the data.

		Returns:
		--------
		list
			Vector indexed by node id containing the node id of mate or networkit.none if unmatched
		"""
		return self._this.getVector()

cdef extern from "<networkit/matching/Matcher.hpp>":

	cdef cppclass _Matcher "NetworKit::Matcher"(_Algorithm):
		_Matcher(const _Graph _G) except +
		_Matching getMatching() except +

cdef class Matcher(Algorithm):
	""" Abstract base class for matching algorithms """
	cdef Graph G

	def __init__(self, *args, **namedargs):
		if type(self) == Matcher:
			raise RuntimeError("Instantiation of abstract base class")

	def getMatching(self):
		"""
		Returns the matching.

		Returns:
		--------
		networkit.Matching
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return Matching().setThis((<_Matcher*>(self._this)).getMatching())


cdef extern from "<networkit/matching/PathGrowingMatcher.hpp>":

	cdef cppclass _PathGrowingMatcher "NetworKit::PathGrowingMatcher"(_Matcher):
		_PathGrowingMatcher(_Graph) except +
		_PathGrowingMatcher(_Graph, vector[double]) except +

cdef class PathGrowingMatcher(Matcher):
	"""
	Path growing matching algorithm as described by  Hougardy and Drake.
	Computes an approximate maximum weight matching with guarantee 1/2.

	Parameters:
	-----------
	G : networkit.Graph
		The input graph (undirected and with no self-loops).
	edgeScores : list, optional
		list with edgeScores.
	"""
	def __cinit__(self, Graph G not None, edgeScores=None):
		self.G = G
		if edgeScores is not None:
			self._this = new _PathGrowingMatcher(G._this, edgeScores)
		else:
			self._this = new _PathGrowingMatcher(G._this)

cdef extern from "<networkit/matching/SuitorMatcher.hpp>":
	cdef cppclass _SuitorMather "NetworKit::SuitorMatcher"(_Matcher):
		_SuitorMather(_Graph, bool_t, bool_t) except +

cdef class SuitorMatcher(Matcher):
	"""
	Computes a 1/2-approximation of the maximum (weighted) matching of an undirected graph using
	the Suitor algorithm from Manne and Halappanavar presented in "New Effective Multithreaded
	Matching Algorithms", IPDPS 2014. The algorithm has two versions: SortSuitor (faster, but
	only works on graphs with adjacency lists sorted by non-increasing edge weight) and Suitor
	(works on generic graphs). If using SortSuitor, call nk.graphtools.sortEdgesByWeight(G, True)
	to sort the adjacency lists by non-increasing edge weight.

	Parameters:
	-----------
	G : networkit.Graph
		The input graph, must be undirected.
	sortSuitor : bool
		If True uses the SortSuitor version, otherwise it uses Suitor.
	checkSortedEdges : bool
		If True and sortSuitor is True it checks whether the adjacency lists
		of the input graph are sorted by non-increasing edge weight. If they are not it raises
		a RuntimeError.
	"""
	def __cinit__(self, Graph G not None, sortSuitor = True, checkSortedEdges = False):
		self.G = G
		self._this = new _SuitorMather(G._this, sortSuitor, checkSortedEdges)
