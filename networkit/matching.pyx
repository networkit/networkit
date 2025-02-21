# distutils: language=c++

from libcpp.vector cimport vector

from .dynbase cimport _DynAlgorithm
from .dynbase import DynAlgorithm
from .dynamics cimport _GraphEvent, GraphEvent
from .helpers import stdstring
from .structures cimport count, node

cdef class Matching:
	""" 
	Matching(z=0)
	
	Implements a graph matching. Create a new matching data structure for `z` elements.

		Parameters
		----------
		z : int, optional
			Maximum number of nodes. Default: 0
	"""
	def __cinit__(self, index z=0):
		self._this = move(_Matching(z))

	cdef setThis(self,  _Matching& other):
		swap[_Matching](self._this,  other)
		return self

	def match(self, node u, node v):
		"""
		match(u, v)
		
		Set two nodes u and v as each others matching partners.

		Parameters
		----------
		u : int
			The first node.
		v : int
			The second node.
		"""
		self._this.match(u,v)

	def unmatch(self, node u,  node v):
		"""
		unmatch(u, v)

		Reset the two nodes u and v to unmatched.

		Parameters
		----------
		u : int
			The first node.
		v : int
			The second node.
		"""
		self._this.unmatch(u, v)

	def isMatched(self, node u):
		"""
		isMatched(u)

		Check if node u is matched.

		Parameters
		----------
		u : int
			The input node.

		Returns
		-------
		bool
			True if u has a matching partner.
		"""
		return self._this.isMatched(u)

	def areMatched(self, node u, node v):
		"""
		areMatched(u, v)

		Check if the nodes u and v are matched to each other.

		Parameters
		----------
		u : int
			The first node.
		v : int
			The second node.

		Returns
		-------
		bool
			True if node u and v are matched to each other.
		"""
		return self._this.areMatched(u,v)

	def isProper(self, Graph G):
		"""
		isProper(G)

		Check whether this is a proper matching in the graph, i.e. there is an edge between
		each pair and if the matching partner of v is u then the matching partner of u is v.

		Parameters
		----------
		G : networkit.Graph
			Graph to be checked.
	
		Returns
		-------
		bool
			True if this is a proper matching.
		"""
		return self._this.isProper(G._this)

	def size(self, Graph G):
		"""
		size(G)

		Get the number of edges in this matching.

		Parameters
		----------
		G : networkit.Graph
			The input graph.

		Returns
		-------
		int
			Total number of edges in the matching.
		"""
		return self._this.size(G._this)

	def mate(self, node v):
		"""
		mate(v)

		Get the matched neighbor of v if it exists, otherwise +inf.

		Parameters
		----------
		v : int
			The input node.

		Returns
		-------
		int
			Matching partner of v if it exists, otherwise +inf.
		"""
		return self._this.mate(v)

	def weight(self, Graph G):
		"""
		weight(G)

		Get total weight of edges in this matching.

		Parameters
		----------
		G : networkit.Graph
			The corresponding graph.

		Returns
		-------
		float
			Total weight of edges in this matching.
		"""
		return self._this.weight(G._this)

	def toPartition(self, Graph G):
		"""
		toPartition(G)

		Convert matching to a Partition object where matched nodes belong to the same subset
		and unmatched nodes belong to a singleton subset.

		Parameters
		----------
		G : networkit.Graph
			The corresponding graph.

		Returns
		-------
		networkit.Partition
			The resulting partition.
		"""
		return Partition().setThis(self._this.toPartition(G._this))

	def getVector(self):
		""" 
		getVector()
		
		Get the vector storing the data.

		Returns
		-------
		list(int or None)
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
		getMatching()

		Returns the matching.

		Returns
		-------
		networkit.Matching
			Current Matching of graph.
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
	PathGrowingMatcher(G, edgeScores)
	
	Path growing matching algorithm as described by  Hougardy and Drake.
	Computes an approximate maximum weight matching with guarantee 1/2.

	Parameters:
	-----------
	G : networkit.Graph
		The input graph (undirected and with no self-loops).
	edgeScores : list(float), optional
		List with edgeScores.
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
	SuitorMatcher(G, sortSuitor=True, checkSortedEdges=False)

	Computes a 1/2-approximation of the maximum (weighted) matching of an undirected graph using
	the Suitor algorithm from Manne and Halappanavar presented in "New Effective Multithreaded
	Matching Algorithms", IPDPS 2014. The algorithm has two versions: SortSuitor (faster, but
	only works on graphs with adjacency lists sorted by non-increasing edge weight) and Suitor
	(works on generic graphs). If using SortSuitor, call nk.graphtools.sortEdgesByWeight(G, True)
	to sort the adjacency lists by non-increasing edge weight.

	Parameters
	----------
	G : networkit.Graph
		The input graph, must be undirected.
	sortSuitor : bool, optional
		If True uses the SortSuitor version, otherwise it uses Suitor. Default: True
	checkSortedEdges : bool, optional
		If True and sortSuitor is True it checks whether the adjacency lists
		of the input graph are sorted by non-increasing edge weight. If they are not it raises
		a RuntimeError. Default: False
	"""
	def __cinit__(self, Graph G not None, sortSuitor = True, checkSortedEdges = False):
		self.G = G
		self._this = new _SuitorMather(G._this, sortSuitor, checkSortedEdges)

cdef class BMatching:
	""" 
	BMatching(b)
	
	Implements a graph matching. Create a new matching data structure for `z` elements.

		Parameters
		----------
		b : list(int)
			List containing the b-values for all nodes.
	"""
	def __cinit__(self, Graph G = None, second = None):
		if G is not None:
			self._G = G
			self._this = move(_BMatching(G._this, <vector[count]> second))
		else:
			self._this = move(_BMatching())

	cdef setThis(self, _BMatching& other):
		swap[_BMatching](self._this,  other)
		return self

	def isProper(self):
		"""
		isProper()

		Check whether this is a proper matching in the graph, i.e. there is an edge between
		each pair and if the matching partner of v is u then the matching partner of u is v.
	
		Returns
		-------
		bool
			True if this is a proper matching.
		"""
		return self._this.isProper()

	def match(self, node u, node v):
		"""
		match(u, v)
		
		Set two nodes u and v as each others matching partners.

		Parameters
		----------
		u : int
			The first node.
		v : int
			The second node.
		"""
		self._this.match(u,v)

	def unmatch(self, node u,  node v):
		"""
		unmatch(u, v)

		Reset the two nodes u and v to unmatched.

		Parameters
		----------
		u : int
			The first node.
		v : int
			The second node.
		"""
		self._this.unmatch(u, v)

	def isUnmatched(self, node u):
		"""
		isMatched(u)

		Check if node u is unmatched.

		Parameters
		----------
		u : int
			The input node.

		Returns
		-------
		bool
			True if u has no matching partner.
		"""
		return self._this.isUnmatched(u)

	def areMatched(self, node u, node v):
		"""
		areMatched(u, v)

		Check if the nodes u and v are matched to each other.

		Parameters
		----------
		u : int
			The first node.
		v : int
			The second node.

		Returns
		-------
		bool
			True if node u and v are matched to each other.
		"""
		return self._this.areMatched(u,v)



	def size(self):
		"""
		size()

		Get the number of edges in this matching.

		Returns
		-------
		int
			Total number of edges in the matching.
		"""
		return self._this.size()

	def weight(self):
		"""
		weight()

		Get total weight of edges in this matching.

		Returns
		-------
		float
			Total weight of edges in this matching.
		"""
		return self._this.weight()

	def getMatches(self):
		"""
		getMatches()

		Get the set of matches for each node.
		
		Returns
		-------
		list(set)
			Returns the set of matches for each node.
		"""
		return self._this.getMatches()

	def getB(self):
		"""
		getB()

		Get the b-value for each node.
		
		Returns
		-------
		list(int)
			A list, containing the b-values (int) for each node.
		"""
		return self._this.getB()

	def reset(self):
		"""
		reset()

		Removes all entries from the b-matching data structure
		"""
		return self._this.reset()

cdef extern from "<networkit/matching/BMatcher.hpp>":

	cdef cppclass _BMatcher "NetworKit::BMatcher"(_Algorithm):
		_BMatcher(const _Graph _G, vector[count] b) except +
		_BMatching getBMatching() except +

cdef class BMatcher(Algorithm):
	""" Abstract base class for matching algorithms """
	cdef Graph _G

	def __init__(self, *args, **namedargs):
		if type(self) == BMatcher:
			raise RuntimeError("Instantiation of abstract base class")

	def getBMatching(self):
		"""
		getMatching()

		Returns the matching.

		Returns
		-------
		networkit.BMatching
			Current b-matching of graph.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return BMatching().setThis((<_BMatcher*>(self._this)).getBMatching())

cdef extern from "<networkit/matching/BSuitorMatcher.hpp>":
	cdef cppclass _BSuitorMatcher "NetworKit::BSuitorMatcher"(_BMatcher):
		_BSuitorMatcher(_Graph, vector[count] b) except +
		_BSuitorMatcher(_Graph, count b) except +
		void buildBMatching() except +

cdef class BSuitorMatcher(BMatcher):
	"""
	BSuitorMatcher(G, second)

    Computes a 1/2-approximate maximum weight b-matching of an undirected weighted Graph @c G
    using the sequential b-Suitor algorithm published by Khan et al. in "Efficient Approximation
    Algorithms For Weighted B-Matching", SIAM Journal on Scientific Computing, Vol. 38, Iss. 5
    (2016).

	Parameters
	----------
	G : networkit.Graph
		The input graph, must be undirected.
    second : The second parameter contains information about b-values of nodes. This can either be given as
             single value (int), indicating the same b-value for all nodes. Other options: list of b-values
             for each node, str-based path to file containing b-values.
	"""

	def __cinit__(self, Graph G, second):

		self._G = G
		if isinstance(second, list):
			self._this = new _BSuitorMatcher(G._this, <vector[count]> second)
		elif isinstance(second, int):
			self._this = new _BSuitorMatcher(G._this, <count> second)
		else:
			raise Exception("Error: the second parameter must be either an int (global b-value), a list of ints (single b-values for all nodes) or a path to the file, containing b-values for every node.")

	def buildBMatching(self):
		"""
		buildBMatching()

    	Creates the b-matching for given graph G. Function run() automatically invokes
    	buildMatching. After invoking buildBMatching(), use getBMatching() to retrieve the resulting
    	b-matching.
		"""
		(<_BSuitorMatcher*>(self._this)).buildBMatching()

cdef extern from "<networkit/matching/DynamicBSuitorMatcher.hpp>":
	cdef cppclass _DynamicBSuitorMatcher "NetworKit::DynamicBSuitorMatcher"(_BSuitorMatcher, _DynAlgorithm):
		_DynamicBSuitorMatcher(_Graph, vector[count] b) except +
		_DynamicBSuitorMatcher(_Graph, count b) except +

cdef class DynamicBSuitorMatcher(BSuitorMatcher, DynAlgorithm):
	""" 
	DynamicBSuitorMatcher(G, second)
	
    Implementation from the algorithm from "A Fully-dynamic Approximation Algorithm for Maximum 
    Weight b-Matchings in Graphs" from Proceedings of The Thirteenth International Conference on 
    Complex Networks and their Applications 2024 by Fabian Brandt-Tumescheit, Frieda Gerharz and
    Henning Meyerhenke. The algorithm dynamically updates the b-matching based on the b-Suitor 
    algorithm by Khan et al. The solution is the same as a complete recomputation of the b-Suitor
    b-matching.

	Parameters
	----------
	G : networkit.Graph
		An unweighted graph.
    second : The second parameter contains information about b-values of nodes. This can either be given as
             single value (int), indicating the same b-value for all nodes. Other options: list of b-values
             for each node, str-based path to file containing b-values.
	"""

	def __cinit__(self, Graph G, second):
		cdef string c_path
		
		self._G = G
		if isinstance(second, list):
			self._this = new _DynamicBSuitorMatcher(G._this, <vector[count]> second)
		elif isinstance(second, int):
			self._this = new _DynamicBSuitorMatcher(G._this, <count> second)
		else:
			raise Exception("Error: the second parameter must be either an int (global b-value), a list of ints (single b-values for all nodes) or a path to the file, containing b-values for every node.")