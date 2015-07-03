""" This module contains algorithms for the sparsification, i.e. edge filtering, of networks. """

__author__ = "Gerd Lindner"

from _NetworKit import ChibaNishizekiTriangleCounter, SimmelianJaccardAttributizer, GlobalThresholdFilter, LocalSimilarityAttributizer, MultiscaleAttributizer, SimmelianOverlapAttributizer, RandomAttributizer, LocalDegreeAttributizer, ForestFireAttributizer, \
	EdgeAttributeAsWeight, EdgeAttributeLinearizer, JaccardSimilarityAttributizer, LocalFilterAttributizer, AdamicAdarDistance, ChanceCorrectedTriangleAttributizer, NodeNormalizedTriangleAttributizer, TriangleCounter, RandomEdgeAttributizer, ChungLuAttributizer, ChibaNishizekiQuadrangleCounter, GeometricMeanAttributizer, \
	EdgeAttributeNormalizer, EdgeAttributeBlender, PrefixJaccardCoefficient, SCANStructuralSimilarityAttributizer

# local imports
from . import community

_ABS_ZERO = 1e-7

class Sparsifier(object):
	""" Abstract base class representing a graph sparsification algorithm that
	uses only one parameter to determine the degree of filtering. """

	def getAttribute(self, G):
		""" Returns an edge attribute. (To be implemented by derived class)

		Keyword arguments:
		G -- the input graph
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
		"""Returns a sparsified version of the given input graph.

		Keyword arguments:
		G -- the input graph
		parameter -- a parameter value that determines the degree of sparsification
		attribute -- (optional) a previously calculated edge attribute. If none is provided, we will try to calculate it.
		"""

		if attribute is None:
			attribute = self.getAttribute(G)

		return self._getSparsifiedGraph(G, parameter, attribute)

	def getSparsifiedGraphOfSize(self, G, edgeRatio, attribute=None):
		"""This is a convenience function that applies an appropriate parameterization
		algorithm (if available) to obtain a parameter value that yields a sparsified
		graph of the desired size and then calls getSparsifiedGraph(...) with that parameter value.

		Keyword arguments:
		G -- the input graph
		edgeRatio -- the target edge ratio
		attribute -- (optional) a previously calculated edge attribute. If none is provided,
		we will try to calculate it.
		"""
		if attribute is None:
			attribute = self.getAttribute(G)

		parameter = self.getParameter(G, edgeRatio, attribute)

		return self.getSparsifiedGraph(G, parameter, attribute)

	def getParameter(self, G, edgeRatio, attribute=None):
		""" This is a convenience function that applies an appropriate parameterization
		algorithm (if available) to obtain a parameter value that yields a sparsified
		graph of the desired size. """
		if attribute is None:
			attribute = self.getAttribute(G)

		paramAlgorithm = self._getParameterizationAlgorithm()
		parameter = paramAlgorithm.parameterize(self, G, attribute, edgeRatio)
		return parameter

class SimpleParameterization:
	""" A parameterization algorithm representds an algorithm that, given a graph
	and a sparsifier, calculates a parameter value such that a desired edge ratio is met.
	The SimpleParameterization strategy simply returns the input edgeRatio as parameterization
	result. """

	def parameterize(self, algorithm, G, attribute, edgeRatio):
		""" Parameterize the given sparsifier for the given input graph with the
		given target edge ratio. (To be implemented by derived class.)

		Keyword arguments:
		algorithm -- the sparsification algorithm
		G -- the input graph
		attribute -- precalculated edge attribute
		edgeRatio -- target edge ratio the resulting parameter value should yield
		"""
		return edgeRatio

class BinarySearchParameterization:
	""" Parameterizes a sparsification algorithm using binary search. """

	def __init__(self, sizeIncreasesWithParameter, lowerParameterBound, upperParameterBound, maxSteps):
		""" Creates a new instance of a binary search parameterizer.

		Keyword arguments:
		sizeIncreasesWithParameter -- set to True if the size of the sparsified graph increases with increasing parameter value
		lowerParameterBound -- lower bound of the parameter domain (inclusive)
		upperParameterBound -- upper bound of the parameter domain (inclusive)
		maxSteps -- the maximum number of steps to perform during binary search
		"""
		self.sizeIncreasesWithParameter = sizeIncreasesWithParameter
		self.lowerParameterBound = lowerParameterBound
		self.upperParameterBound = upperParameterBound
		self.maxSteps = maxSteps

	def parameterize(self, algorithm, G, attribute, edgeRatio):
		lowerBound = self.lowerParameterBound
		upperBound = self.upperParameterBound
		estimation = self.lowerParameterBound
		bestParameter = self.lowerParameterBound
		minDistance = self.upperParameterBound

		for i in range(0, self.maxSteps):
			estimation = (lowerBound + upperBound) / 2.0
			backbone = algorithm._getSparsifiedGraph(G, estimation, attribute)
			currentEdgeRatio = backbone.numberOfEdges() / G.numberOfEdges()

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
	""" Parameterizes a sparsification algorithm using complete search
	(applicable only to algorithms which take as input a parameter from a small
	set of possible values) """

	def __init__(self, lowerParameterBound, upperParameterBound):
		""" Creates a new instance of a complete search parameterizer.

		Keyword arguments:
		lowerParameterBound -- lower bound of the parameter domain (inclusive, integer)
		upperParameterBound -- upper bound of the parameter domain (inclusive, integer) """
		self.lowerParameterBound = lowerParameterBound
		self.upperParameterBound = upperParameterBound

	def parameterize(self, algorithm, G, attribute, edgeRatio):
		bestParameter = self.lowerParameterBound
		bestRatio = 0.0
		minDistance = 100.0

		for i in range(self.lowerParameterBound, self.upperParameterBound + 1):
			backbone = algorithm._getSparsifiedGraph(G, i, attribute)
			currentEdgeRatio = backbone.numberOfEdges() / G.numberOfEdges()

			distance = abs(currentEdgeRatio - edgeRatio)
			if distance < minDistance and abs(currentEdgeRatio) > _ABS_ZERO:
				minDistance = distance
				bestParameter = i
				bestRatio = currentEdgeRatio

		return bestParameter

def getRankAttribute(attribute, reverse = False):
	""" Takes as input an attribute (node or edge) and returns an attribute where
	each node is assigned its rank among all others according to the attribute values.
	The node/edge with lowest input value is assigned 0, the one with second-lowest
	value 1, and so on.

	Keyword arguments:
	attribute -- the input node/edge attribute
	reverse -- reverses the ranking, if set to True

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

class SimmelianBackboneParametric(Sparsifier):

	""" An implementation of the Parametric variant of the Simmelian Backbones
	 introduced by Nick et al. """

	def __init__(self, maxRank):
		self.maxRank = maxRank

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		chiba = ChibaNishizekiTriangleCounter(G)
		triangles = chiba.getAttribute()
		so = SimmelianOverlapAttributizer(G, triangles, self.maxRank)
		a_so = so.getAttribute()
		return a_so

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return CompleteSearchParameterization(0, self.maxRank)

class SimmelianBackboneNonParametric(Sparsifier):

	""" An implementation of the Non-parametric variant of the Simmelian Backbones
	introduced by Nick et al. """

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum jaccard filter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		chiba = ChibaNishizekiTriangleCounter(G)
		triangles = chiba.getAttribute()
		sj = SimmelianJaccardAttributizer(G, triangles)
		a_sj = sj.getAttribute()
		return a_sj

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class SimmelianMultiscaleBackbone(Sparsifier):

	""" Multiscale Backbone that uses triangle counts as input edge weight. """

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the maximum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		chiba = ChibaNishizekiTriangleCounter(G)
		triangles = chiba.getAttribute()
		ms = MultiscaleAttributizer(G, triangles)
		a_ms = ms.getAttribute()
		return a_ms

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class DegreeMultiscaleBackbone(Sparsifier):
	""" Multiscale Backbone that uses node degrees (mapped to edges) as input edge weight. """

	def __init__(self, degsToAttrValue):
		"""
		Creates a new instance of the Degree Multiscale sparsifier.
		Keyword arguments:
		attrType -- For each edge (x,y) the following edge value will be used.
			0: max(d(x), d(y))
			1: min(d(x), d(y))
			2: avg(d(x), d(y))
		"""
		self.degsToAttrValue = degsToAttrValue

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the maximum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		inputAttribute = [0] * G.upperEdgeIdBound()
		for (x,y) in G.edges():
			inputAttribute[G.edgeId(x,y)] = self.degsToAttrValue(G.degree(x), G.degree(y))

		ms = MultiscaleAttributizer(G, inputAttribute)
		a_ms = ms.getAttribute()
		return a_ms

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class LocalSimilarityBackbone(Sparsifier):

	""" An implementation of the Local Similarity sparsification approach introduced by Satuluri et al. """

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		chiba = ChibaNishizekiTriangleCounter(G)
		triangles = chiba.getAttribute()
		attributizer = LocalSimilarityAttributizer(G, triangles)
		a_ls = attributizer.getAttribute()
		return a_ls

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class MultiscaleBackbone(Sparsifier):

	""" An implementation of the Multiscale backbone approach introduced by Serrano et al. """

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		inputAttribute = [0.0] * G.upperEdgeIdBound()
		for edge in G.edges():
			edgeId = G.edgeId(edge[0], edge[1])
			inputAttribute[edgeId] = G.weight(edge[0], edge[1])

		attributizer = MultiscaleAttributizer(G, inputAttribute)
		attribute = attributizer.getAttribute()
		return attribute

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class RandomBackbone(Sparsifier):

	""" Random Edge sampling. Edges to keep in the backbone are selected uniformly at random. """

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		attributizer = RandomAttributizer(G)
		a_r = attributizer.getAttribute()
		return a_r

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, False)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return SimpleParameterization()

class ForestFireBackbone(Sparsifier):

	""" A variant of the Forest Fire sparsification approach proposed by Leskovec et al. """

	def __init__(self, burnProbability, targetBurntRatio):
		""" Creates a new instance of the Edge Forest Fire backbone algorithm.

		Keyword arguments:
		burnProbability -- the probability that the neighbor of a burnt node gets burnt as well.
		edgeRatio -- the fires will stop when a edgeRatio * edgeCount edges have been burnt.
		"""
		self.burnProbability = burnProbability
		self.targetBurntRatio = targetBurntRatio

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the maximum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		attributizer = ForestFireAttributizer(G, self.burnProbability, self.targetBurntRatio)
		return attributizer.getAttribute()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		 return BinarySearchParameterization(False, 0.0, 1.0, 20)

class LocalDegreeBackbone(Sparsifier):

	""" An implementation of the Local Degree backbone algorithm. """

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		attributizer_ld = LocalDegreeAttributizer(G)
		a_ld = attributizer_ld.getAttribute()
		return a_ld

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class TriangleBackbone(Sparsifier):
	"""  Allows for global filtering with respect to triangle counts. """

	def getAttribute(self, G):
		""" Returns the triangle counts of the input graph. """
		attributizer_t = TriangleCounter(G)
		a_t = attributizer_t.getAttribute()
		return a_t

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		raise NotImplementedError("parameterization method not yet implemented.")

class ModularityPartitionAttributizer():

	"""  """

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		cdAlgo = community.PLM(G, refine=True, turbo=True)
		cdAlgo.run()
		partition = cdAlgo.getPartition()

		def together(u, v):
			if (partition[u] == partition[v]):
				return 1.0
			else:
				return 0.0

		# FIXME: with respect to performance, the following is wrong on so many levels - don't try this at home
		edgeScores = [None for i in range(G.upperEdgeIdBound())]
		for (u, v) in G.edges():
			edgeScores[G.edgeId(u, v)] = together(u, v)
		return edgeScores

class ConstantAttributizer():
	""" Assigns as an attribute the same value to each edge (for sanity checks) """

	def __init__(self, constValue = 1.0):
		""" Creates a new instance of an attributizer that always
		 returns the given value as edge attribute value.
		"""
		self.constValue = constValue

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the constant value given
		in the constructor.
		"""
		return self.constValue
