""" This module contains algorithms for the sparsification, i.e. edge filtering, of networks. The methods are described in detail in the publication
	` Structure-Preserving Sparsification of Social Networks` at http://arxiv.org/abs/1505.00564
"""

__author__ = "Gerd Lindner"

from _NetworKit import ChibaNishizekiTriangleEdgeScore, GlobalThresholdFilter, LocalSimilarityScore, MultiscaleScore, SimmelianOverlapScore, RandomEdgeScore, LocalDegreeScore, ForestFireScore, \
	EdgeScoreAsWeight, EdgeScoreLinearizer, LocalFilterScore, AdamicAdarDistance, ChanceCorrectedTriangleScore, TriangleEdgeScore, RandomNodeEdgeScore, ChibaNishizekiQuadrangleEdgeScore, GeometricMeanScore, \
	EdgeScoreNormalizer, EdgeScoreBlender, PrefixJaccardScore, SCANStructuralSimilarityScore, JaccardSimilarityAttributizer

# local imports
from . import community
from . import distance

_ABS_ZERO = 1e-7

class Sparsifier(object):
	""" Abstract base class representing a graph sparsification algorithm that
	uses only one parameter to determine the degree of filtering. """

	def scores(self, G):
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
			attribute = self.scores(G)

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
			attribute = self.scores(G)

		parameter = self.getParameter(G, edgeRatio, attribute)

		return self.getSparsifiedGraph(G, parameter, attribute)

	def getParameter(self, G, edgeRatio, attribute=None):
		""" This is a convenience function that applies an appropriate parameterization
		algorithm (if available) to obtain a parameter value that yields a sparsified
		graph of the desired size. """
		if attribute is None:
			attribute = self.scores(G)

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
			sparsified = algorithm._getSparsifiedGraph(G, i, attribute)
			currentEdgeRatio = sparsified.numberOfEdges() / G.numberOfEdges()

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

class SimmelianSparsifierParametric(Sparsifier):

	""" An implementation of the Parametric variant of the Simmelian Sparsifiers
	 introduced by Nick et al. """

	def __init__(self, maxRank):
		self.maxRank = maxRank

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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

	""" An implementation of the Non-parametric variant of the Simmelian Sparsifiers
	introduced by Nick et al. """

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the minimum jaccard filter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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
	""" An implementation of the Simmelian Sparsifiers based on quadrangles. """

	def scores(self, G):
		"""
		Returns an edge scoring attribute that can be used for global filtering.

		Keyword arguments:
		G -- the input graph
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

	""" Multiscale Sparsifier that uses triangle counts as input edge weight. """

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the maximum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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
	""" Multiscale Sparsifier that uses node degrees (mapped to edges) as input edge weight. """

	def __init__(self, degsToAttrValue):
		"""
		Creates a new instance of the Degree Multiscale sparsifier.
		Keyword arguments:
		degsToAttrValue -- function that maps two node degrees to an edge score.
		"""
		self.degsToAttrValue = degsToAttrValue

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the maximum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		inputAttribute = [0] * G.upperEdgeIdBound()
		for (x,y) in G.edges():
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
	""" An implementation of the Jaccard Similarity sparsification approach introduced by Satuluri et al. """

	def scores(self, G):
		""" Returns the jaccard coefficient of the neighborhoods of the two incident nodes """
		triangles = TriangleEdgeScore(G).run().scores()
		return JaccardSimilarityAttributizer(G, triangles).getAttribute()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)


class LocalSimilaritySparsifier(Sparsifier):

	""" An implementation of the Local Similarity sparsification approach introduced by Satuluri et al. """

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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

	""" An implementation of the Multiscale backbone approach introduced by Serrano et al. """

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
		"""

		inputAttribute = [0.0] * G.upperEdgeIdBound()
		for edge in G.edges():
			edgeId = G.edgeId(edge[0], edge[1])
			inputAttribute[edgeId] = G.weight(edge[0], edge[1])

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

	""" Random Edge sampling. Edges to keep in the sparsified graph are selected uniformly at random. """

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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
	""" [TODO not yet documented] """

	def __init__(self, above = True):
		self.above = above

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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

	""" A variant of the Forest Fire sparsification approach proposed by Leskovec et al. """

	def __init__(self, burnProbability, targetBurntRatio):
		""" Creates a new instance of the Edge Forest Fire sparsification algorithm.

		Keyword arguments:
		burnProbability -- the probability that the neighbor of a burnt node gets burnt as well.
		edgeRatio -- the fires will stop when a edgeRatio * edgeCount edges have been burnt.
		"""
		self.burnProbability = burnProbability
		self.targetBurntRatio = targetBurntRatio

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the maximum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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

	""" An implementation of the Local Degree sparsification algorithm. """

	def scores(self, G):
		""" Returns an edge score that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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

	""" A sparsifiier dervived from 'SCAN: a structural clustering algorithm for networks' """

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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
	"""  Allows for global filtering with respect to triangle counts. """

	def scores(self, G):
		""" Returns the triangle counts of the input graph. """
		triangleScore = TriangleEdgeScore(G)
		triangleScore.run()
		return triangleScore.scores()

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		raise NotImplementedError("parameterization method not yet implemented.")

class AlgebraicDistanceSparsifier(Sparsifier):
	""" Allows for global filtering with respect to (inverted) algebraic distances. """

	def __init__(self, numberSystems=10, numberIterations=30, omega=0.5, norm=0):
		self.numberSystems = numberSystems
		self.numberIterations = numberIterations
		self.omega = omega
		self.norm = norm

	def scores(self, G):
		""" Returns the inverted algebraic distance score of the input graph. """
		algDist = distance.AlgebraicDistance(G, self.numberSystems, self.numberIterations, self.omega, self.norm, withEdgeScores=True)
		algDist.preprocess()
		return [1.0 - d for d in algDist.getEdgeScores()]

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(G, attribute, parameter, True)
		return gf.calculate()

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

class LocalSparsifier(Sparsifier):
       def __init__(self, sparsifier):
               self.sparsifier = sparsifier

       def scores(self, G):
               """ Returns an edge attribute that holds for each edge 1 - the minimum parameter value
               such that the edge is contained in the sparsified graph.

               Note that - like for all sparsifiers - edges with the highest score are the most important ones.

               Keyword arguments:
               G -- the input graph
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

	"""  """

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.

		Keyword arguments:
		G -- the input graph
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
	""" Assigns as an attribute the same value to each edge (for sanity checks) """

	def __init__(self, constValue = 1.0):
		""" Creates a new instance of an attributizer that always
		 returns the given value as edge attribute value.
		"""
		self.constValue = constValue

	def scores(self, G):
		""" Returns an edge attribute that holds for each edge the constant value given
		in the constructor.
		"""
		return self.constValue
