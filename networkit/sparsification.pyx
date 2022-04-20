# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .distance import AdamicAdarDistance, JaccardSimilarityAttributizer
from .structures cimport count, node, index, edgeid, edgeweight
from . import community
from . import distance


cdef extern from "<algorithm>" namespace "std":
	vector[double] move(vector[double])
	vector[bool_t] move(vector[bool_t])

cdef extern from "<networkit/edgescores/EdgeScore.hpp>":

	cdef cppclass _EdgeScore "NetworKit::EdgeScore"[T](_Algorithm):
		_EdgeScore(const _Graph& G) except +
		vector[T] scores() except +
		T score(edgeid eid) except +
		T score(node u, node v) except +

cdef class EdgeScore(Algorithm):
	"""
	Abstract base class for an edge score.
	"""
	cdef Graph _G

	cdef bool_t isDoubleValue(self):
		raise RuntimeError("Implement in subclass")

	def __init__(self, *args, **namedargs):
		if type(self) == EdgeScore:
			raise RuntimeError("Error, you may not use EdgeScore directly, use a sub-class instead")

	def __dealloc__(self):
		self._G = None # just to be sure the graph is deleted

	def score(self, u, v = None):
		"""
		score(u, v=None)

		Get a vector containing the score for a certain edge connecting node u and v.

		Parameters
		----------
		u : int
			The first node.
		v : int, optional
			The second node. Default: None

		Returns
		-------
		list(int) or list(float)
			List with scores. Scores are either int-values (unweighted) or float-values (weighted).
		"""
		if v is None:
			if self.isDoubleValue():
				return (<_EdgeScore[double]*>(self._this)).score(u)
			else:
				return (<_EdgeScore[count]*>(self._this)).score(u)
		else:
			if self.isDoubleValue():
				return (<_EdgeScore[double]*>(self._this)).score(u, v)
			else:
				return (<_EdgeScore[count]*>(self._this)).score(u, v)

	def scores(self):
		"""
		scores()

		Get a vector containing the score for each edge in the graph.

		Returns
		-------
		list(int) or list(float)
			List with scores. Scores are either int-values (unweighted) or float-values (weighted).
		"""
		if self.isDoubleValue():
			return (<_EdgeScore[double]*>(self._this)).scores()
		else:
			return (<_EdgeScore[count]*>(self._this)).scores()


cdef extern from "<networkit/edgescores/ChibaNishizekiTriangleEdgeScore.hpp>":

	cdef cppclass _ChibaNishizekiTriangleEdgeScore "NetworKit::ChibaNishizekiTriangleEdgeScore"(_EdgeScore[count]):
		_ChibaNishizekiTriangleEdgeScore(const _Graph& G) except +

cdef class ChibaNishizekiTriangleEdgeScore(EdgeScore):
	"""
	ChibaNishizekiTriangleEdgeScore(G)

	Calculates for each edge the number of triangles it is embedded in.

	Parameters
	----------
	G : networkit.Graph
		The graph to count triangles on.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _ChibaNishizekiTriangleEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return False

cdef extern from "<networkit/edgescores/ChibaNishizekiQuadrangleEdgeScore.hpp>":

	cdef cppclass _ChibaNishizekiQuadrangleEdgeScore "NetworKit::ChibaNishizekiQuadrangleEdgeScore"(_EdgeScore[count]):
		_ChibaNishizekiQuadrangleEdgeScore(const _Graph& G) except +

cdef class ChibaNishizekiQuadrangleEdgeScore(EdgeScore):
	"""
	ChibaNishizekiQuadrangleEdgeScore(G)

	Calculates for each edge the number of quadrangles (circles of length 4) it is embedded in.

	Parameters
	----------
	G : networkit.Graph
		The graph to count quadrangles on.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _ChibaNishizekiQuadrangleEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return False

cdef extern from "<networkit/edgescores/TriangleEdgeScore.hpp>":

	cdef cppclass _TriangleEdgeScore "NetworKit::TriangleEdgeScore"(_EdgeScore[double]):
		_TriangleEdgeScore(const _Graph& G) except +

cdef class TriangleEdgeScore(EdgeScore):
	"""
	TriangleEdgeScore(G)

	Triangle counting.

	Parameters
	----------
	G : networkit.Graph
		The graph to count triangles on.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _TriangleEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return False

cdef extern from "<networkit/edgescores/EdgeScoreLinearizer.hpp>":

	cdef cppclass _EdgeScoreLinearizer "NetworKit::EdgeScoreLinearizer"(_EdgeScore[double]):
		_EdgeScoreLinearizer(const _Graph& G, const vector[double]& attribute, bool_t inverse) except +

cdef class EdgeScoreLinearizer(EdgeScore):
	"""
	EdgeScoreLinearizer(G, a, inverse=False)

	Linearizes a score such that values are evenly distributed between 0 and 1.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	a : list(float)
		Edge score that shall be linearized.
	inverse : bool, optional
		Set to True in order to inverse the resulting score. Default: False
	"""
	cdef vector[double] _score

	def __cinit__(self, Graph G, vector[double] score, inverse = False):
		self._G = G
		self._score = score
		self._this = new _EdgeScoreLinearizer(G._this, self._score, inverse)

	cdef bool_t isDoubleValue(self):
		return True


cdef extern from "<networkit/edgescores/EdgeScoreNormalizer.hpp>":

	cdef cppclass _EdgeScoreNormalizer "NetworKit::EdgeScoreNormalizer"[T](_EdgeScore[double]):
		_EdgeScoreNormalizer(const _Graph&, vector[T]&, bool_t inverse, double lower, double upper) except +

cdef class EdgeScoreNormalizer(EdgeScore):
	"""
	EdgeScoreNormalizer(G, score, inverse=False, lower=0.0, uppper=1.0)

	Normalize an edge score such that it is in a certain range.

	Parameters
	----------
	G : networkit.Graph
		The graph the edge score is defined on.
	score : list(float)
		The edge score to normalize.
	inverse: bool, optional
		Set to True in order to inverse the resulting score. Default: False
	lower : float, optional
		Lower bound of the target range. Default: 0.0
	upper : float, optional
		Upper bound of the target range. Default: 0.0
	"""
	cdef vector[double] _inScoreDouble
	cdef vector[count] _inScoreCount

	def __cinit__(self, Graph G not None, score, bool_t inverse = False, double lower = 0.0, double upper = 1.0):
		self._G = G
		try:
			self._inScoreDouble = <vector[double]?>score
			self._this = new _EdgeScoreNormalizer[double](G._this, self._inScoreDouble, inverse, lower, upper)
		except TypeError:
			try:
				self._inScoreCount = <vector[count]?>score
				self._this = new _EdgeScoreNormalizer[count](G._this, self._inScoreCount, inverse, lower, upper)
			except TypeError:
				raise TypeError("score must be either a vector of integer or float")

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/edgescores/EdgeScoreBlender.hpp>":

	cdef cppclass _EdgeScoreBlender "NetworKit::EdgeScoreBlender"(_EdgeScore[double]):
		_EdgeScoreBlender(const _Graph&, const vector[double]&, const vector[double]&, const vector[bool_t]&) except +

cdef class EdgeScoreBlender(EdgeScore):
	"""
	EdgeScoreBlender(G, attribute0, attribute1, selection)

	Blends two attribute vectors, the value is chosen depending on the supplied bool vector

	Parameters
	----------
	G : networkit.Graph
		The graph for which the attribute shall be blended
	attribute0 : list(float)
		The first attribute (chosen for selection[eid] == False)
	attribute1 : list(float)
		The second attribute (chosen for selection[eid] == True)
	selection : list(bool)
		The selection vector.
	"""
	cdef vector[double] _attribute0
	cdef vector[double] _attribute1
	cdef vector[bool_t] _selection

	def __cinit__(self, Graph G not None, vector[double] attribute0, vector[double] attribute1, vector[bool_t] selection):
		self._G = G
		self._attribute0 = move(attribute0)
		self._attribute1 = move(attribute1)
		self._selection = move(selection)

		self._this = new _EdgeScoreBlender(G._this, self._attribute0, self._attribute1, self._selection)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/edgescores/GeometricMeanScore.hpp>":

	cdef cppclass _GeometricMeanScore "NetworKit::GeometricMeanScore"(_EdgeScore[double]):
		_GeometricMeanScore(const _Graph& G, const vector[double]& a) except +

cdef class GeometricMeanScore(EdgeScore):
	"""
	GeometricMeanScore(G, a)

	Normalizes the given edge attribute by the geometric average of the sum of the attributes of the incident edges of the incident nodes.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	a : list(float)
		Edge attribute that shall be normalized.
	"""
	cdef vector[double] _attribute

	def __cinit__(self, Graph G, vector[double] attribute):
		self._G = G
		self._attribute = attribute
		self._this = new _GeometricMeanScore(G._this, self._attribute)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/edgescores/EdgeScoreAsWeight.hpp>":

	cdef cppclass _EdgeScoreAsWeight "NetworKit::EdgeScoreAsWeight":
		_EdgeScoreAsWeight(const _Graph& G, const vector[double]& score, bool_t squared, edgeweight offset, edgeweight factor) except +
		_Graph calculate() except +

cdef class EdgeScoreAsWeight:
	"""
	EdgeScoreAsWeight(G, score, square, offset, factor)

	Assigns an edge score as edge weight of a graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to assign edge weights to.
	score : list(float)
		The input edge score.
	squared : bool
		Edge weights will be squared if set to True.
	offset : float
		This offset will be added to each edge weight.
	factor : float
		Each edge weight will be multiplied by this factor.
	"""

	cdef _EdgeScoreAsWeight* _this
	cdef Graph _G
	cdef vector[double] _score

	def __cinit__(self, Graph G, vector[double] score, bool_t squared, edgeweight offset, edgeweight factor):
		self._G = G
		self._score = score
		self._this = new _EdgeScoreAsWeight(G._this, self._score, squared, offset, factor)

	def __dealloc__(self):
		if self._this is not NULL:
			del self._this
			self._this = NULL

	def getWeightedGraph(self):
		"""
		getWeightedGraph()

		Returns
		-------
		networkit.Graph
			The weighted result graph.
		"""
		return Graph(0).setThis(self._this.calculate())

cdef extern from "<networkit/sparsification/SimmelianOverlapScore.hpp>":

	cdef cppclass _SimmelianOverlapScore "NetworKit::SimmelianOverlapScore"(_EdgeScore[double]):
		_SimmelianOverlapScore(const _Graph& G, const vector[count]& triangles, count maxRank) except +

cdef class SimmelianOverlapScore(EdgeScore):
	"""
	SimmelianOverlapScore(G, triangles, maxRank)

	An implementation of the parametric variant of Simmelian Backbones. Calculates
	for each edge the minimum parameter value such that the edge is still contained in
	the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Simmelian Backbone algorithm to.
	triangles : list(int)
		Previously calculated edge triangle counts on G.
	maxRank: int
		Maximum rank that is considered for overlap calculation.
	"""
	cdef vector[count] _triangles

	def __cinit__(self, Graph G, vector[count] triangles, count maxRank):
		self._G = G
		self._triangles = triangles
		self._this = new _SimmelianOverlapScore(G._this, self._triangles, maxRank)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/edgescores/PrefixJaccardScore.hpp>":

	cdef cppclass _PrefixJaccardScore "NetworKit::PrefixJaccardScore<double>"(_EdgeScore[double]):
		_PrefixJaccardScore(const _Graph& G, const vector[double]& a) except +

cdef class PrefixJaccardScore(EdgeScore):
	cdef vector[double] _attribute

	def __cinit__(self, Graph G, vector[double] attribute):
		self._G = G
		self._attribute = attribute
		self._this = new _PrefixJaccardScore(G._this, self._attribute)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/MultiscaleScore.hpp>":

	cdef cppclass _MultiscaleScore "NetworKit::MultiscaleScore"(_EdgeScore[double]):
		_MultiscaleScore(const _Graph& G, const vector[double]& a) except +

cdef class MultiscaleScore(EdgeScore):
	"""
	MultiscaleScore(G, attribute)

	An implementation of the Multiscale Backbone. Calculates for each edge the minimum
	parameter value such that the edge is still contained in the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Multiscale algorithm to.
	attribute : list(float)
		The edge attribute the Multiscale algorithm is to be applied to.
	"""

	cdef vector[double] _attribute

	def __cinit__(self, Graph G, vector[double] attribute):
		self._G = G
		self._attribute = attribute
		self._this = new _MultiscaleScore(G._this, self._attribute)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/RandomEdgeScore.hpp>":

	cdef cppclass _RandomEdgeScore "NetworKit::RandomEdgeScore"(_EdgeScore[double]):
		_RandomEdgeScore(const _Graph& G) except +

cdef class RandomEdgeScore(EdgeScore):
	"""
	RandomEdgeScore(G)

	Generates a random edge attribute. Each edge is assigned a random value in [0,1].

	Parameters
	----------
	G : networkit.Graph
		The graph to calculate the Random Edge attribute for.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _RandomEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/LocalSimilarityScore.hpp>":

	cdef cppclass _LocalSimilarityScore "NetworKit::LocalSimilarityScore"(_EdgeScore[double]):
		_LocalSimilarityScore(const _Graph& G, const vector[count]& triangles) except +

cdef class LocalSimilarityScore(EdgeScore):
	"""
	LocalSimilarityScore(G, triangles)

	An implementation of the Local Simlarity sparsification approach.
	This attributizer calculates for each edge the maximum parameter value
	such that the edge is still contained in the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Local Similarity algorithm to.
	triangles : list(int)
		Previously calculated edge triangle counts.
	"""
	cdef vector[count] _triangles

	def __cinit__(self, Graph G, vector[count] triangles):
		self._G = G
		self._triangles = triangles
		self._this = new _LocalSimilarityScore(G._this, self._triangles)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/ForestFireScore.hpp>":

	cdef cppclass _ForestFireScore "NetworKit::ForestFireScore"(_EdgeScore[double]):
		_ForestFireScore(const _Graph& G, double pf, double tebr) except +

cdef class ForestFireScore(EdgeScore):
	"""
	ForestFireScore(G, pf, tebr)

	A variant of the Forest Fire sparsification approach that is based on random walks.
	This attributizer calculates for each edge the minimum parameter value
	such that the edge is still contained in the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Forest Fire algorithm to.
	pf : float
		The probability for neighbor nodes to get burned aswell.
	tebr : float
		The Forest Fire will burn until tebr * numberOfEdges edges have been burnt.
	"""

	def __cinit__(self, Graph G, double pf, double tebr):
		self._G = G
		self._this = new _ForestFireScore(G._this, pf, tebr)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/LocalDegreeScore.hpp>":

	cdef cppclass _LocalDegreeScore "NetworKit::LocalDegreeScore"(_EdgeScore[double]):
		_LocalDegreeScore(const _Graph& G) except +

cdef class LocalDegreeScore(EdgeScore):
	"""
	LocalDegreeScore(G)

	The LocalDegree sparsification approach is based on the idea of hub nodes.
	This attributizer calculates for each edge the maximum parameter value
	such that the edge is still contained in the sparsified graph.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Local Degree  algorithm to.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _LocalDegreeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/RandomNodeEdgeScore.hpp>":

	cdef cppclass _RandomNodeEdgeScore "NetworKit::RandomNodeEdgeScore"(_EdgeScore[double]):
		_RandomNodeEdgeScore(const _Graph& G) except +

cdef class RandomNodeEdgeScore(EdgeScore):
	"""
	RandomNodeEdgeScore(G)

	Random Edge sampling. This attributizer returns edge attributes where
	each value is selected uniformly at random from [0,1].

	Parameters
	----------
	G : networkit.Graph
		The graph to calculate the Random Edge attribute for.
	"""
	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _RandomNodeEdgeScore(G._this)

	cdef bool_t isDoubleValue(self):
		return True

ctypedef fused DoubleInt:
	int
	double

cdef extern from "<networkit/sparsification/LocalFilterScore.hpp>":

	cdef cppclass _LocalFilterScoreDouble "NetworKit::LocalFilterScore<double>"(_EdgeScore[double]):
		_LocalFilterScoreDouble(const _Graph& G, const vector[double]& a, bool_t logarithmic) except +

	cdef cppclass _LocalFilterScoreInt "NetworKit::LocalFilterScore<int>"(_EdgeScore[count]):
		_LocalFilterScoreInt(const _Graph& G, const vector[double]& a, bool_t logarithmic) except +

cdef class LocalFilterScore(EdgeScore):
	"""
	LocalFilterScore(G, a, logarithmic=True)

	Local filtering edge scoring. Edges with high score are more important.

	Edges are ranked locally, the top d^e (logarithmic, default) or 1+e*(d-1) edges (non-logarithmic) are kept.
	For equal attribute values, neighbors of low degree are preferred.
	If bothRequired is set (default: false), both neighbors need to indicate that they want to keep the edge.

	Parameters
	----------
	G : networkit.Graph
		The input graph
	a : list(float)
		The input attribute according to which the edges shall be fitlered locally.
	logarithmic : bool, optional
		If the score shall be logarithmic in the rank (then d^e edges are kept), linear otherwise. Default: True
	"""
	cdef vector[double] _a

	def __cinit__(self, Graph G, vector[double] a, bool_t logarithmic = True):
		self._G = G
		self._a = a
		self._this = new _LocalFilterScoreDouble(G._this, self._a, logarithmic)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/ChanceCorrectedTriangleScore.hpp>":

	cdef cppclass _ChanceCorrectedTriangleScore "NetworKit::ChanceCorrectedTriangleScore"(_EdgeScore[double]):
		_ChanceCorrectedTriangleScore(const _Graph& G, const vector[count]& triangles) except +

cdef class ChanceCorrectedTriangleScore(EdgeScore):
	"""
	ChanceCorrectedTriangleScore(G, triangles)

	Divide the number of triangles per edge by the expected number of triangles given a random edge distribution.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	triangles : list(int)
		Triangle count.
	"""
	cdef vector[count] _triangles

	def __cinit__(self, Graph G, vector[count] triangles):
		self._G = G
		self._triangles = triangles
		self._this = new _ChanceCorrectedTriangleScore(G._this, self._triangles)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/SCANStructuralSimilarityScore.hpp>":

	cdef cppclass _SCANStructuralSimilarityScore "NetworKit::SCANStructuralSimilarityScore"(_EdgeScore[double]):
		_SCANStructuralSimilarityScore(_Graph G, const vector[count]& triangles) except +

cdef class SCANStructuralSimilarityScore(EdgeScore):
	"""
	SCANStructuralSimilarityScore(G, triangles)

	An implementation of the SCANStructuralSimilarityScore algorithm.

	Parameters
	----------
	G : networkit.Graph
		The graph to apply the Local Similarity algorithm to.
	triangles : list(int)
		Previously calculated edge triangle counts.
	"""
	cdef vector[count] _triangles

	def __cinit__(self, Graph G, vector[count] triangles):
		self._G = G
		self._triangles = triangles
		self._this = new _SCANStructuralSimilarityScore(G._this, self._triangles)

	cdef bool_t isDoubleValue(self):
		return True

cdef extern from "<networkit/sparsification/GlobalThresholdFilter.hpp>":

	cdef cppclass _GlobalThresholdFilter "NetworKit::GlobalThresholdFilter":
		_GlobalThresholdFilter(const _Graph& G, const vector[double]& a, double alpha, bool_t above) except +
		_Graph calculate() except +

cdef class GlobalThresholdFilter:
	"""
	GlobalThresholdFilter(G, attribute, e, above)

	Calculates a sparsified graph by filtering globally using a constant threshold value
	and a given edge attribute.

	Parameters
	----------
	G : networkit.Graph
		The graph to sparsify.
	attribute : list(float)
		The edge attribute to consider for filtering.
	e : float
		Threshold value.
	above : bool
		If set to True (False), all edges with an attribute value equal to or above (below)
		will be kept in the sparsified graph.
	"""
	cdef _GlobalThresholdFilter* _this
	cdef Graph _G
	cdef vector[double] _attribute

	def __cinit__(self, Graph G not None, vector[double] attribute, double e, bool_t above):
		self._G = G
		self._attribute = attribute
		self._this = new _GlobalThresholdFilter(G._this, self._attribute, e, above)

	def __dealloc__(self):
		del self._this

	def calculate(self):
		return Graph().setThis(self._this.calculate())

_ABS_ZERO = 1e-7

class Sparsifier(object):
	""" Abstract base class representing a graph sparsification algorithm that
	uses only one parameter to determine the degree of filtering. """

	def scores(self, G):
		"""
		scores(G)

		Returns an edge attribute.

		Parameters
		----------
		G : networkit.Graph
			The input graph.

		Returns
		-------
		list(float)
			Edge scores.
		"""
		raise NotImplementedError

	def _getSparsifiedGraph(self, G, parameter, attribute):
		""" The actual implementation of the sparsification.
		(To be implemented in the derived class.)

		Keyword arguments:
		G -- the input graph
		edgeRatio -- the target edge ratio
		attribute -- a previously calculated edge attribute.
		"""
		raise NotImplementedError

	def _getParameterizationAlgorithm(self):
		""" Returns an appropriate parameterization algorithm for this sparsifier.
		(To be implemented in the derived class.)
		"""
		return SimpleParameterization()

	def getSparsifiedGraph(self, G, parameter, attribute=None):
		"""
		getSparsifiedGraph(sG, parameter, attribute=None)
		
		Returns a sparsified version of the given input graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		parameter : float
			A parameter value that determines the degree of sparsification
		attribute : list(float), optional
			A previously calculated edge attribute. If none is provided, we will try to calculate it. Default: None
		"""

		if attribute is None:
			attribute = self.scores(G)

		return self._getSparsifiedGraph(G, parameter, attribute)

	def getSparsifiedGraphOfSize(self, G, edgeRatio, attribute=None):
		"""
		getSparsifiedGraphOfSize(G, edgeRatio, attribute=None)
		
		This is a convenience function that applies an appropriate parameterization
		algorithm (if available) to obtain a parameter value that yields a sparsified
		graph of the desired size and then calls getSparsifiedGraph(...) with that parameter value.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		edgeRatio : float
			The target edge ratio.
		attribute : list(float), optional
			A previously calculated edge attribute. If none is provided, we will try to calculate it. Default: None

		Returns
		-------
		networkit.Graph
			The sparsified graph.
		"""
		if attribute is None:
			attribute = self.scores(G)

		parameter = self.getParameter(G, edgeRatio, attribute)

		return self.getSparsifiedGraph(G, parameter, attribute)

	def getParameter(self, G, edgeRatio, attribute=None):
		""" 
		getParameter(G, edgeRatio, attribute=None)
		
		This is a convenience function that applies an appropriate parameterization
		algorithm (if available) to obtain a parameter value that yields a sparsified
		graph of the desired size.
		
		Parameters
		----------
		G : networkit.Graph
			The input graph.
		edgeRatio : float
			The target edge ratio.
		attribute : list(float), optional
			A previously calculated edge attribute. If none is provided, we will try to calculate it. Default: None
		"""
		if attribute is None:
			attribute = self.scores(G)

		paramAlgorithm = self._getParameterizationAlgorithm()
		parameter = paramAlgorithm.parameterize(self, G, attribute, edgeRatio)
		return parameter

class SimpleParameterization:
	""" 
	SimpleParameterization()
	
	A parameterization algorithm representds an algorithm that, given a graph
	and a sparsifier, calculates a parameter value such that a desired edge ratio is met.
	The SimpleParameterization strategy simply returns the input edgeRatio as parameterization
	result. 
	"""

	def parameterize(self, algorithm, G, attribute, edgeRatio):
		""" 
		parameterize(algorithm, G, attribute, edgeRatio)
		
		Parameterize the given sparsifier for the given input graph with the
		given target edge ratio. (To be implemented by derived class.)

		Parameters
		----------
		algorithm : networkit.sparsification.Sparsifier
			The sparsification algorithm.
		G : networkit.Graph
			The input graph.
		attribute : list(float)
			An input edge attribute.
		edgeRatio : float
			The target edge ratio.
		"""
		return edgeRatio

class BinarySearchParameterization:
	""" 
	BinarySearchParameterization(sizeIncreasesWithParameter, lowerParameterBound, upperParameterBound, maxSteps)
	
	Parameterizes a sparsification algorithm using binary search.

	Parameters
	----------
	sizeIncreasesWithParameter : bool
		Set to True if the size of the sparsified graph increases with increasing parameter value.
	lowerParameterBound : float
		Lower bound of the parameter domain (inclusive).
	upperParameterBound : float
		Upper bound of the parameter domain (inclusive).
	maxSteps : int
		The maximum number of steps to perform during binary search.
	"""

	def __init__(self, sizeIncreasesWithParameter, lowerParameterBound, upperParameterBound, maxSteps):

		self.sizeIncreasesWithParameter = sizeIncreasesWithParameter
		self.lowerParameterBound = lowerParameterBound
		self.upperParameterBound = upperParameterBound
		self.maxSteps = maxSteps

	def parameterize(self, algorithm, G, attribute, edgeRatio):
		""" 
		parameterize(algorithm, G, attribute, edgeRatio)
		
		Parameterize the given sparsifier for the given input graph with the
		given target edge ratio. (To be implemented by derived class.)

		Parameters
		----------
		algorithm : networkit.sparsification.Sparsifier
			The sparsification algorithm.
		G : networkit.Graph
			The input graph.
		attribute : list(float)
			An input edge attribute.
		edgeRatio : float
			The target edge ratio.
		"""
		lowerBound = self.lowerParameterBound
		upperBound = self.upperParameterBound
		estimation = self.lowerParameterBound
		bestParameter = self.lowerParameterBound
		minDistance = self.upperParameterBound

		for i in range(0, self.maxSteps):
			estimation = (lowerBound + upperBound) / 2.0
			sparsified = algorithm._getSparsifiedGraph(G, estimation, attribute)
			currentEdgeRatio = sparsified.numberOfEdges() / G.numberOfEdges()

			distance = abs(currentEdgeRatio - edgeRatio)

			if distance < minDistance and abs(currentEdgeRatio) > _ABS_ZERO:
				minDistance = distance
				bestParameter = estimation

				#"Exact" hit?
				if abs(currentEdgeRatio - edgeRatio) < _ABS_ZERO:
					break;

			increase = (currentEdgeRatio < edgeRatio)
			if not self.sizeIncreasesWithParameter:
				increase = not increase

			if increase:
				lowerBound = estimation
			else:
				upperBound = estimation

		return bestParameter

class CompleteSearchParameterization:
	""" 
	CompleteSearchParameterization(lowerParameterBound, upperParameterBound)
	
	Parameterizes a sparsification algorithm using complete search
	(applicable only to algorithms which take as input a parameter from a small
	set of possible values).
	
	Parameters
	----------
	lowerParameterBound : int
		Lower bound of the parameter domain (inclusive).
	upperParameterBound : int
		Upper bound of the parameter domain (inclusive).
	"""

	def __init__(self, lowerParameterBound, upperParameterBound):
		self.lowerParameterBound = lowerParameterBound
		self.upperParameterBound = upperParameterBound

	def parameterize(self, algorithm, G, attribute, edgeRatio):
		""" 
		parameterize(algorithm, G, attribute, edgeRatio)
		
		Parameterize the given sparsifier for the given input graph with the
		given target edge ratio. (To be implemented by derived class.)

		Parameters
		----------
		algorithm : networkit.sparsification.Sparsifier
			The sparsification algorithm.
		G : networkit.Graph
			The input graph.
		attribute : list(float)
			An input edge attribute.
		edgeRatio : float
			The target edge ratio.
		"""
		bestParameter = self.lowerParameterBound
		bestRatio = 0.0
		minDistance = 100.0

		for i in range(self.lowerParameterBound, self.upperParameterBound + 1):
			sparsified = algorithm._getSparsifiedGraph(G, i, attribute)
			currentEdgeRatio = sparsified.numberOfEdges() / G.numberOfEdges()

			distance = abs(currentEdgeRatio - edgeRatio)
			if distance < minDistance and abs(currentEdgeRatio) > _ABS_ZERO:
				minDistance = distance
				bestParameter = i
				bestRatio = currentEdgeRatio

		return bestParameter

def getRankAttribute(attribute, reverse = False):
	""" 
	getRankAttribute(attribute, reverse=False)
	
	Takes as input an attribute (node or edge) and returns an attribute where
	each node is assigned its rank among all others according to the attribute values.
	The node/edge with lowest input value is assigned 0, the one with second-lowest
	value 1, and so on.

	Parameters
	----------
	attribute : list(float)
		The input node/edge attribute
	reverse : bool, optional
		Reverses the ranking, if set to True. Default: False
	"""

	#Example input: [0.1, 0.05, 0.9, 0.2], ascending
	#Example output: [1, 0, 3, 2]

	_attribute = zip([x for x in range(0, len(attribute))], attribute)
	_attribute = sorted(_attribute, key=lambda x: x[1], reverse=reverse)

	_index = 0
	result = [0] * len(attribute)
	for (i, v) in _attribute:
		result[i] = _index
		_index = _index + 1

	return result

class SimmelianSparsifierParametric(Sparsifier):
	""" 
	SimmelianSparsifierParametric()

	An implementation of the Parametric variant of the Simmelian Sparsifiers
	introduced by Nick et al.
	"""

	def __init__(self, maxRank):
		self.maxRank = maxRank

	def scores(self, G):
		""" 
		scores(G)		

		Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		triangles = TriangleEdgeScore(G).run().scores()

		simmelianOverlap = SimmelianOverlapScore(G, triangles, self.maxRank)
		simmelianOverlap.run()
		return simmelianOverlap.scores()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return CompleteSearchParameterization(0, self.maxRank)

class SimmelianSparsifierNonParametric(Sparsifier):
	"""
	SimmelianSparsifierNonParametric()
	
	An implementation of the Non-parametric variant of the Simmelian Sparsifiers
	introduced by Nick et al. 
	"""

	def scores(self, G):
		""" 
		scores(G)		

		Returns an edge attribute that holds for each edge the minimum jaccard filter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		triangles = TriangleEdgeScore(G).run().scores()
		a_sj = PrefixJaccardScore(G, triangles).run().scores()

		return a_sj

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class QuadrilateralSimmelianSparsifier(Sparsifier):
	""" 
	QuadrilateralSimmelianSparsifier()	

	An implementation of the Simmelian Sparsifiers based on quadrangles.
	"""

	def scores(self, G):
		"""
		scores(G)

		Returns an edge scoring attribute that can be used for global filtering.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		quadrangles = ChibaNishizekiQuadrangleEdgeScore(G).run().scores()
		meanQuadrangles = GeometricMeanScore(G, quadrangles).run().scores()
		quadranglePrefixJaccard = PrefixJaccardScore(G, meanQuadrangles).run().scores()
		return quadranglePrefixJaccard

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class SimmelianMultiscaleSparsifier(Sparsifier):
	"""
	SimmelianMultiscaleSparsifier()	

	Multiscale Sparsifier that uses triangle counts as input edge weight.
	"""

	def scores(self, G):
		""" 
		scores(G)		

		Returns an edge attribute that holds for each edge the maximum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		triangles = TriangleEdgeScore(G).run().scores()
		ms = MultiscaleScore(G, triangles)
		ms.run()
		a_ms = ms.scores()
		return a_ms

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class DegreeMultiscaleSparsifier(Sparsifier):
	""" 
	DegreeMultiscaleSparsifier(degsToAttrValue)

	Multiscale Sparsifier that uses node degrees (mapped to edges) as input edge weight.

	Parameters
	----------
	degsToAttrValue : function
		Function that maps two node degrees to an edge score.
	"""

	def __init__(self, degsToAttrValue):
		self.degsToAttrValue = degsToAttrValue

	def scores(self, G):
		""" 
		scores(G)		

		Returns an edge attribute that holds for each edge the maximum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""

		inputAttribute = [0] * G.upperEdgeIdBound()
		for x, y in G.iterEdges():
			inputAttribute[G.edgeId(x,y)] = self.degsToAttrValue(G.degree(x), G.degree(y))

		ms = MultiscaleScore(G, inputAttribute)
		ms.run()
		a_ms = ms.scores()
		return a_ms

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class JaccardSimilaritySparsifier(Sparsifier):
	""" 
	JaccardSimilaritySparsifier()

	An implementation of the Jaccard Similarity sparsification approach introduced by Satuluri et al.
	"""

	def scores(self, G):
		""" 
		scores(G)		

		Returns the jaccard coefficient of the neighborhoods of the two incident nodes.
		
		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		triangles = TriangleEdgeScore(G).run().scores()
		return JaccardSimilarityAttributizer(G, triangles).getAttribute()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)


class LocalSimilaritySparsifier(Sparsifier):
	"""
	LocalSimilaritySparsifier()

	An implementation of the Local Similarity sparsification approach introduced by Satuluri et al.
	"""

	def scores(self, G):
		"""
		scores(G)

		Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		triangles = TriangleEdgeScore(G).run().scores()
		localSimScore = LocalSimilarityScore(G, triangles)
		localSimScore.run()
		return localSimScore.scores()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class MultiscaleSparsifier(Sparsifier):
	""" 
	MultiscaleSparsifier()

	An implementation of the Multiscale backbone approach introduced by Serrano et al.
	"""

	def scores(self, G):
		""" 
		scores(G)

		Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""

		inputAttribute = [0.0] * G.upperEdgeIdBound()
		for u, v, w in G.iterEdgesWeights():
			edgeId = G.edgeId(u, v)
			inputAttribute[edgeId] = w

		scorer = MultiscaleScore(G, inputAttribute)
		scorer.run()
		scores = scorer.scores()
		return scores

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class RandomEdgeSparsifier(Sparsifier):
	"""
	RandomEdgeSparsifier()

	Random Edge sampling. Edges to keep in the sparsified graph are selected uniformly at random.
	"""

	def scores(self, G):
		""" 
		scores(G)

		Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""

		reScore = RandomEdgeScore(G)
		reScore.run()
		return reScore.scores()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, False)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return SimpleParameterization()

class RandomNodeEdgeSparsifier(Sparsifier):
	"""
	RandomNodeEdgeSparsifier()

	Random Edge sampling. Edges and nodes to keep in the sparsified graph are selected uniformly at random.

	Parameters
	----------
	above : bool, optional
		If set to True (False), all nodes and edges with an attribute value equal to or above (below)
		will be kept in the sparsified graph. Default: True
	"""

	def __init__(self, above = True):
		self.above = above

	def scores(self, G):
		""" 
		scores(G)

		Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""

		rneScore = RandomNodeEdgeScore(G)
		rneScore.run()
		return rneScore.scores()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, self.above)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization((not self.above), 0.0, 1.0, 20)

class ForestFireSparsifier(Sparsifier):
	"""
	ForestFireSparsifier(burnProbability, targetBurntRatio)
	
	A variant of the Forest Fire sparsification approach proposed by Leskovec et al.
	
	Parameters
	----------
	burnProbability : float
		The probability that the neighbor of a burnt node gets burnt as well.
	edgeRatio : float
		The fires will stop when a edgeRatio * edgeCount edges have been burnt.
	"""

	def __init__(self, burnProbability, targetBurntRatio):
		self.burnProbability = burnProbability
		self.targetBurntRatio = targetBurntRatio

	def scores(self, G):
		""" 
		scores(G)
		
		Returns an edge attribute that holds for each edge the maximum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""

		ffScore = ForestFireScore(G, self.burnProbability, self.targetBurntRatio)
		ffScore.run()
		return ffScore.scores()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class LocalDegreeSparsifier(Sparsifier):
	""" 
	LocalDegreeSparsifier()
	
	An implementation of the Local Degree sparsification algorithm.
	"""

	def scores(self, G):
		""" 
		scores(G)
		
		Returns an edge score that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""

		localDegree = LocalDegreeScore(G)
		localDegree.run()
		localDegreeScore = localDegree.scores()
		return localDegreeScore

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class SCANSparsifier(Sparsifier):
	""" 
	SCANSparsifier()
	
	A sparsifiier dervived from 'SCAN: a structural clustering algorithm for networks'
	"""

	def scores(self, G):
		"""
		scores(G)

		Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		a_triangles = TriangleEdgeScore(G).run().scores()

		scanScore = SCANStructuralSimilarityScore(G, a_triangles)
		scanScore.run()

		return scanScore.scores()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class TriangleSparsifier(Sparsifier):
	""" 
	TriangleSparsifier()

	Allows for global filtering with respect to triangle counts.
	"""

	def scores(self, G):
		""" 
		scores(G)

		Returns the triangle counts of the input graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		triangleScore = TriangleEdgeScore(G)
		triangleScore.run()
		return triangleScore.scores()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		raise NotImplementedError("parameterization method not yet implemented.")

class AlgebraicDistanceSparsifier(Sparsifier):
	""" 
	AlgebraicDistanceSparsifier(numberSystems=10, numberIterations=30, omega=0.5, norm=0)
	
	Allows for global filtering with respect to (inverted) algebraic distances.
	
	Parameters
	----------
	numberSystems : int, optional
	 	Number of vectors/systems used for algebraic iteration. Default: 10
	numberIterations : int, optional
	 	Number of iterations in each system. Default: 30
	omega : float, optional
	 	Attenuation factor in [0,1] influencing convergence speed. Default: 0.5
	norm : int, optional
		The norm factor of the extended algebraic distance. Default: 0
	"""

	def __init__(self, numberSystems=10, numberIterations=30, omega=0.5, norm=0):
		self.numberSystems = numberSystems
		self.numberIterations = numberIterations
		self.omega = omega
		self.norm = norm

	def scores(self, G):
		""" 
		scores(G)

		Returns the inverted algebraic distance score of the input graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		algDist = distance.AlgebraicDistance(G, self.numberSystems, self.numberIterations, self.omega, self.norm, withEdgeScores=True)
		algDist.preprocess()
		return [1.0 - d for d in algDist.getEdgeScores()]

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class LocalSparsifier(Sparsifier):
	"""
	LocalSparsifier(sparsifier)

	Sparsification algorithm using an input sparsifier and passing the result of
	this algorithm to networkit.sparsification.LocalFilterScore for choosing
	edges, which should be present in the sparsified graph.
	"""
	def __init__(self, sparsifier):
		self.sparsifier = sparsifier

	def scores(self, G):
		""" 
		scores(G)		

		Returns an edge attribute that holds for each edge 1 - the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		originalScores = self.sparsifier.scores(G)
		localFilterScore = LocalFilterScore(G, originalScores)
		localFilterScore.run()

		return localFilterScore.scores()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)


class ModularityPartitionScore():
	"""
	ModularityPartitionScore()

	Sparsification algorithm using networkit.community.PLM algorithm to detect
	communities and select only edges, where both ends are member of the same
	community.
	"""

	def scores(self, G):
		""" 
		scores(G)

		Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		cdAlgo = community.PLM(G, par="none randomized", refine=True, turbo=True)
		cdAlgo.run()
		partition = cdAlgo.getPartition()

		def together(u, v):
			if (partition[u] == partition[v]):
				return 1.0
			else:
				return 0.0

		edgeScores = [None for i in range(G.upperEdgeIdBound())]
		G.forEdges(lambda u, v, w, eid: edgeScores.__setitem__(eid, together(u, v)))
		return edgeScores

class ConstantScore():
	""" 
	ConstantScore(constValue=1.0)

	Assigns as an attribute the same value to each edge (for sanity checks).

	Parameters
	----------
	constValue : float, optional
		Value to assign to each edge. Default: 1.0
	"""

	def __init__(self, constValue = 1.0):
		self.constValue = constValue

	def scores(self, G):
		"""
		scores(G)

		Returns an edge attribute that holds for each edge the constant value given
		in the constructor.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		return self.constValue
