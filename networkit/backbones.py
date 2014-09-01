""" This module contains algorithms for the sparsification, i.e. edge filtering, of networks. """

__author__ = "Gerd Lindner"

from _NetworKit import ChibaNishizekiTriangleCounter, SimmelianJaccardAttributizer, GlobalThresholdFilter, LocalSimilarityAttributizer, MultiscaleAttributizer, SimmelianOverlapAttributizer, RandomAttributizer, LocalDegreeAttributizer, ForestFireAttributizer, \
	AttributeAsWeight, LinearizeAttribute, JaccardSimilarityAttributizer, LocalFilterAttributizer, TriangleCounter

#SimmelianBackboneParametric, SimmelianBackboneNonParametric, MultiscaleBackbone, LocalSimilarityBackbone, SimmelianMultiscaleBackbone, RandomBackbone

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

	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.
		
		Keyword arguments:
		G -- the input graph
		"""
		
		chiba = ChibaNishizekiTriangleCounter()
		triangles = chiba.getAttribute(G)
		so = SimmelianOverlapAttributizer(10)
		a_so = so.getAttribute(G, triangles)
		return a_so

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, True)
		return gf.calculate(G, attribute)

	def _getParameterizationAlgorithm(self):
		return CompleteSearchParameterization(0, 10)
	
class SimmelianBackboneNonParametric(Sparsifier):
	
	""" An implementation of the Non-parametric variant of the Simmelian Backbones 
	introduced by Nick et al. """
	
	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum jaccard filter value
		such that the edge is contained in the sparsified graph.
		
		Keyword arguments:
		G -- the input graph
		"""
		
		chiba = ChibaNishizekiTriangleCounter()
		triangles = chiba.getAttribute(G)
		sj = SimmelianJaccardAttributizer()
		a_sj = sj.getAttribute(G, triangles)
		return a_sj
	
	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, True)
		return gf.calculate(G, attribute)

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
		
		chiba = ChibaNishizekiTriangleCounter()
		triangles = chiba.getAttribute(G)
		ms = MultiscaleAttributizer()
		a_ms = ms.getAttribute(G, triangles)
		return a_ms
       
	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)
	
	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)

class LocalSimilarityBackbone(Sparsifier):
	
	""" An implementation of the Local Similarity sparsification approach introduced by Satuluri et al. """
	
	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.
		
		Keyword arguments:
		G -- the input graph
		"""
		
		attributizer = LocalSimilarityAttributizer()
		a_ls = attributizer.getAttribute(G, [])
		return a_ls

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)

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
		
		attributizer = MultiscaleAttributizer()
		attribute = attributizer.getAttribute(G, inputAttribute)
		return attribute

	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)
	
	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)
	
class RandomBackbone(Sparsifier):
	
	""" Random Edge sampling. Edges to keep in the backbone are selected uniformly at random. """
	
	def getAttribute(self, G):
		""" Returns an edge attribute that holds for each edge the minimum parameter value
		such that the edge is contained in the sparsified graph.
		
		Keyword arguments:
		G -- the input graph
		"""
		
		attributizer = RandomAttributizer(1.0)
		a_r = attributizer.getAttribute(G, [1.0] * G.upperEdgeIdBound())
		return a_r
		
	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)
	
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
		
		attributizer = ForestFireAttributizer(self.burnProbability, self.targetBurntRatio)
		return attributizer.getAttribute(G, [])
		
	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, True)
		return gf.calculate(G, attribute)
	
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
		
		attributizer_ld = LocalDegreeAttributizer()
		a_ld = attributizer_ld.getAttribute(G, [])
		return a_ld
       
	def _getSparsifiedGraph(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)
	
	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)

