from libc.stdint cimport uint64_t
from libcpp cimport bool as bool_t
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string

ctypedef uint64_t index
ctypedef index node

from .graph cimport _Graph, Graph
from .base cimport _Algorithm, Algorithm
from .structures cimport _Cover, Cover

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes

cdef extern from "<networkit/scd/PageRankNibble.hpp>":

	cdef cppclass _PageRankNibble "NetworKit::PageRankNibble":
		_PageRankNibble(_Graph G, double alpha, double epsilon) except +
		map[node, set[node]] run(set[node] seeds) except +

cdef class PageRankNibble:
	"""
	Produces a cut around a given seed node using the PageRank-Nibble algorithm.
	see Andersen, Chung, Lang: Local Graph Partitioning using PageRank Vectors

	Parameters:
	-----------
	G : networkit.Graph in which the cut is to be produced, must be unweighted.
	alpha : Loop probability of random walk; smaller values tend to produce larger communities.
	epsilon: Tolerance threshold for approximation of PageRank vectors
	"""
	cdef _PageRankNibble *_this
	cdef Graph _G

	def __cinit__(self, Graph G, double alpha, double epsilon):
		self._G = G
		self._this = new _PageRankNibble(G._this, alpha, epsilon)

	def run(self, set[node] seeds):
		"""
		Produces a cut around a given seed node.

		Parameters:
		-----------
		seeds : the seed node ids.
		"""
		return self._this.run(seeds)

cdef extern from "<networkit/scd/GCE.hpp>":

	cdef cppclass _GCE "NetworKit::GCE":
		_GCE(_Graph G, string quality) except +
		map[node, set[node]] run(set[node] seeds) except +

cdef class GCE:
	"""
	Produces a cut around a given seed node using the GCE algorithm.

	Parameters:
	-----------
	G : networkit.Graph in which the cut is to be produced, must be unweighted.
	"""
	cdef _GCE *_this
	cdef Graph _G

	def __cinit__(self, Graph G, quality):
		self._G = G
		self._this = new _GCE(G._this, stdstring(quality))

	def run(self, set[node] seeds):
		"""
		Produces a cut around a given seed node.

		Parameters:
		-----------
		seeds : the seed node ids.
		"""
		return self._this.run(seeds)

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
