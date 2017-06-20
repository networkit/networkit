from libc.stdint cimport uint64_t
from libcpp cimport bool as bool_t
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string

ctypedef uint64_t index
ctypedef index node

from cython.operator import dereference

from .graph cimport _Graph, Graph
from .base cimport _Algorithm, Algorithm
from .structures cimport _Cover, Cover

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes

cdef extern from "<networkit/scd/SelectiveCommunityDetector.hpp>":

	cdef cppclass _SelectiveCommunityDetector "NetworKit::SelectiveCommunityDetector":
		_SelectiveCommunityDetector(_Graph G) except +
		map[node, set[node]] run(set[node] seeds) except +
		set[node] expandOneCommunity(node seed) except +
		set[node] expandOneCommunity(set[node] seeds) except +

cdef class SelectiveCommunityDetector:
	cdef _SelectiveCommunityDetector *_this
	cdef Graph _G

	def __init__(self, *args, **namedargs):
		if type(self) == SelectiveCommunityDetector:
			raise RuntimeError("Error, you may not use SelectiveCommunityDetector directly, use a sub-class instead")

	def __cinit__(self, *args, **namedargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

	def run(self, set[node] seeds):
		"""
		Detect one community for each of the given seed nodes.

		The default implementation calls expandOneCommunity() for each of the seeds.

		Parameters:
		-----------
		seeds : list of nodes
			The list of seeds for which communities shall be detected.

		Returns
		-------
		A dict mapping from seed node to community (as a set of nodes).
		"""
		return self._this.run(seeds)

	def expandOneCommunity(self, seeds):
		"""
		Detect a community for the given seed node(s).

		It expands either a single seed or a whole set of seed nodes into a single community.
		This is useful if you know multiple nodes that should be part of the
		community. This method may throw an exception if the particular algorithm
		does not support multiple seeds but you specified more than one seed.

		Parameters
		----------
		seed : node or set
			The seed(s) to find the community for.

		Returns
		-------
		The found community as a set of nodes.
		"""
		try:
			return self._this.expandOneCommunity(<set[node]?>seeds)
		except TypeError:
			return self._this.expandOneCommunity(<node?>seeds)


cdef extern from "<networkit/scd/PageRankNibble.hpp>":

	cdef cppclass _PageRankNibble "NetworKit::PageRankNibble"(_SelectiveCommunityDetector):
		_PageRankNibble(_Graph G, double alpha, double epsilon) except +

cdef class PageRankNibble(SelectiveCommunityDetector):
	"""
	Produces a cut around a given seed node using the PageRank-Nibble algorithm.
	see Andersen, Chung, Lang: Local Graph Partitioning using PageRank Vectors

	Parameters:
	-----------
	G : networkit.Graph in which the cut is to be produced, must be unweighted.
	alpha : Loop probability of random walk; smaller values tend to produce larger communities.
	epsilon: Tolerance threshold for approximation of PageRank vectors
	"""
	def __cinit__(self, Graph G, double alpha, double epsilon):
		self._G = G
		self._this = new _PageRankNibble(G._this, alpha, epsilon)

cdef extern from "<networkit/scd/GCE.hpp>":

	cdef cppclass _GCE "NetworKit::GCE"(_SelectiveCommunityDetector):
		_GCE(_Graph G, string quality) except +

cdef class GCE(SelectiveCommunityDetector):
	"""
	Produces a cut around a given seed node using the GCE algorithm.
	It greedily adds nodes from the shell to improve community quality.

	Parameters:
	-----------
	G : graph in which the cut is to be produced, must be unweighted.
	Q : string
		The quality function. Supported values: "M" or "L".
	"""
	def __cinit__(self, Graph G, quality):
		self._G = G
		self._this = new _GCE(G._this, stdstring(quality))

cdef extern from "<networkit/scd/CliqueDetect.hpp>":
	cdef cppclass _CliqueDetect "NetworKit::CliqueDetect"(_SelectiveCommunityDetector):
		_CliqueDetect(_Graph G) except +

cdef class CliqueDetect(SelectiveCommunityDetector):
	"""
	The CliqueDetect algorithm. It finds the largest clique in the
	seed node's neighborhood.

	The algorithm can handle weighted graphs. There, the clique
	with the highest sum of internal edge weights is
	returned. This sum includes edge weights to the seed node(s)
	to ensure that cliques that are well-connected to the seed
	node(s) are preferred.

	For multiple seed nodes, the resulting community is a clique
	iff the seeds form a clique. Otherwise, only the added nodes
	form a clique that is fully connected to the seed nodes.

	See also: Hamann, M.; Röhrs, E.; Wagner, D.
	Local Community Detection Based on Small Cliques.
	Algorithms 2017, 10, 90. https://doi.org/10.3390/a10030090

	Parameters:
	-----------
	G : graph in which communities shall be detected.
	"""
	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _CliqueDetect(G._this)

cdef extern from "<networkit/scd/LFMLocal.hpp>":

	cdef cppclass _LFMLocal "NetworKit::LFMLocal"(_SelectiveCommunityDetector):
		_LFMLocal(_Graph G, double alpha) except +

cdef class LFMLocal(SelectiveCommunityDetector):
	"""
	Local version of the LFM algorithm

	This is the local community expansion as introduced in:

	Lancichinetti, A., Fortunato, S., & Kertész, J. (2009).
	Detecting the overlapping and hierarchical community structure in complex networks.
	New Journal of Physics, 11(3), 033015.
	https://doi.org/10.1088/1367-2630/11/3/033015

	Their algorithm detects overlapping communities by repeatedly
	executing this algorithm for a random seed node that has not yet
	been assigned to any community.

	The algorithm has a resolution parameter alpha. A natural choice
	for alpha is 1, the paper states that values below 0.5 usually
	give a community containing the whole graph while values larger
	than 2 recover the smallest communities.

	Parameters:
	-----------
	G : Graph
		graph in which the community shall be found.
	alpha : float
		The resolution parameter
	"""
	def __cinit__(self, Graph G, double alpha = 1.0):
		self._G = G
		self._this = new _LFMLocal(G._this, alpha)

cdef extern from "<networkit/scd/CombinedSCD.hpp>":
	cdef cppclass _CombinedSCD "NetworKit::CombinedSCD"(_SelectiveCommunityDetector):
		_CombinedSCD(_Graph G, _SelectiveCommunityDetector first, _SelectiveCommunityDetector second) except +

cdef class CombinedSCD(SelectiveCommunityDetector):
	cdef SelectiveCommunityDetector _first
	cdef SelectiveCommunityDetector _second

	def __cinit__(self, Graph G, SelectiveCommunityDetector first not None, SelectiveCommunityDetector second not None):
		self._G = G
		self._first = first
		self._second = second
		self._this = new _CombinedSCD(G._this, dereference(first._this), dereference(second._this))

cdef extern from "<networkit/scd/SCDGroundTruthComparison.hpp>":
	cdef cppclass _SCDGroundTruthComparison "NetworKit::SCDGroundTruthComparison"(_Algorithm):
		_SCDGroundTruthComparison(_Graph G, _Cover groundTruth, map[node, set[node]] found, bool_t ignoreSeeds) except +
		map[index, double] getIndividualJaccard() except +
		map[index, double] getIndividualPrecision() except +
		map[index, double] getIndividualRecall() except +
		map[index, double] getIndividualF1() except +

		double getAverageJaccard() except +
		double getAverageF1() except +
		double getAveragePrecision() except +
		double getAverageRecall() except +

cdef class SCDGroundTruthComparison(Algorithm):
	"""
	This class evaluates a set found communities against a ground truth cover. Each found community
	is compared against the communities of the seed node in the ground truth cover.

	For each score, the ground truth community is chosen as comparison that maximizes the score.
	If seeds are not ignored (a parameter of the constructor), then only ground truth communities
	that contain the given seed are used to compare against.

	The calculated scores are:

	Precision: the size of the intersection of found and ground truth community divided by the
	size of the found community, i.e., how much of the found community was an actual match.

	Recall: the size of the intersection of found and ground truth community divided by the
	size of the ground truth community, i.e., how much of the ground truth community was found.

	F1 score: the harmonic mean of precision and recall.

	Jaccard index: the size of the intersection of found and ground truth community divided by the
	size of the union of found and ground truth community.

	For each score, the range of values is between 0 and 1, where 0 is the worst and 1 the best score.

	Parameters:
	-----------
	G : Graph
		The graph to compare on
	groundTruth : Cover
		The ground truth cover
	found : dict
		The found communities, where the keys are the seed nodes and the values are the found nodes.
	ignoreSeeds : bool
		If the seed values shall be ignored, i.e. if any ground truth community is a match
	"""
	cdef Graph _G
	cdef Cover _groundTruth
	cdef map[node, set[node]] _found

	def __cinit__(self, Graph G not None, Cover groundTruth not None, map[node, set[node]] found, bool_t ignoreSeeds):
		self._found = found
		self._G = G
		self._groundTruth = groundTruth
		self._this = new _SCDGroundTruthComparison(G._this, groundTruth._this, self._found, ignoreSeeds)

	def getIndividualJaccard(self):
		"""
		Get the Jaccard index of every found community.

		Returns
		-------
		A dict between seed node and the jaccard index of the seed's community.
		"""
		return (<_SCDGroundTruthComparison*>(self._this)).getIndividualJaccard()
	def getIndividualPrecision(self):
		"""
		Get the precision of every found community.

		Returns
		-------
		A dict between seed node and the precision of the seed's community.
		"""
		return (<_SCDGroundTruthComparison*>(self._this)).getIndividualPrecision()
	def getIndividualRecall(self):
		"""
		Get the recall of every found community.

		Returns
		-------
		A dict between seed node and the recall of the seed's community.
		"""
		return (<_SCDGroundTruthComparison*>(self._this)).getIndividualRecall()
	def getIndividualF1(self):
		"""
		Get the F1 score of every found community.

		Returns
		-------
		A dict between seed node and the F1 score of the seed's community.
		"""
		return (<_SCDGroundTruthComparison*>(self._this)).getIndividualF1()
	def getAverageJaccard(self):
		"""
		Get the (unweighted) average of the jaccard indices of every found community.
		"""
		return (<_SCDGroundTruthComparison*>(self._this)).getAverageJaccard()
	def getAverageF1(self):
		"""
		Get the (unweighted) average of the F1 score of every found community.
		"""
		return (<_SCDGroundTruthComparison*>(self._this)).getAverageF1()
	def getAveragePrecision(self):
		"""
		Get the (unweighted) average of the precision of every found community.
		"""
		return (<_SCDGroundTruthComparison*>(self._this)).getAveragePrecision()
	def getAverageRecall(self):
		"""
		Get the (unweighted) average of the recall of every found community.
		"""
		return (<_SCDGroundTruthComparison*>(self._this)).getAverageRecall()

cdef extern from "<networkit/scd/SetConductance.hpp>":

	cdef cppclass _SetConductance "NetworKit::SetConductance"(_Algorithm):
		_SetConductance(_Graph G, set[node]) except +
		double getConductance() except +

cdef class SetConductance(Algorithm):
	"""
	This class calculates the conductance of a set of nodes, i.e., the weight of all edges
	between the set and the rest of the graph divided by the minimum
	of the volume (the sum of the weighted degrees) of the community and the rest
	of the graph.

	Parameters:
	-----------
	G : Graph
		The graph to calculate the conductance on
	community : set
		The set of nodes to calculate the conductance of
	"""

	cdef Graph _G
	cdef set[node] _community

	def __cinit__(self, Graph G not None, set[node] community):
		self._community = community
		self._G = G
		self._this = new _SetConductance(G._this, self._community)

	def getConductance(self):
		"""
		Get the calculated conductance score.

		Returns
		-------
		The conductance.
		"""
		return (<_SetConductance*>(self._this)).getConductance()

cdef extern from "<networkit/scd/RandomBFS.hpp>":

	cdef cppclass _RandomBFS "NetworKit::RandomBFS"(_SelectiveCommunityDetector):
		_RandomBFS(_Graph G, _Cover C) except +

cdef class RandomBFS(SelectiveCommunityDetector):
	"""
	The random BFS community detection baseline:
	finds a community around a seed node with a given size using a prefix of a
	BFS order that is selected at random among all possible prefixes.

	Parameters:
	-----------
	G : Graph
		Graph in which the community shall be found
	C : Cover
		Ground truth communities to get size information from
	"""
	cdef Cover _C

	def __cinit__(self, Graph G, Cover C):
		self._G = G
		self._C = C
		self._this = new _RandomBFS(G._this, C._this)
