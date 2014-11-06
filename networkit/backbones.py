""" This module contains algorithms for the sparsification of networks, i.e. Backbone algorithms. """

__author__ = "Gerd Lindner"

import abc as _abc
from _NetworKit import ChibaNishizekiTriangleCounter, SimmelianJaccardAttributizer, GlobalThresholdFilter, LocalSimilarityAttributizer, MultiscaleAttributizer, SimmelianOverlapAttributizer, RandomAttributizer, LocalDegreeAttributizer, ForestFireAttributizer

#SimmelianBackboneParametric, SimmelianBackboneNonParametric, MultiscaleBackbone, LocalSimilarityBackbone, SimmelianMultiscaleBackbone, RandomBackbone

_ABS_ZERO = 1e-7

""" Abstract base class representing a graph sparsification algorithm. 
 Sparsification algorithms can be split up into a precalculation step
 (calculation of an attribute) and a global filtering step. """
class BackboneAlgorithm(object):
	__metaclass__ = _abc.ABCMeta
	
	@_abc.abstractmethod
	def getAttribute(self, G):
		"""Returns an edge attribute"""
		pass

	def getBackbone(self, G, parameter, attribute=None):
		"""Returns a sparsified version of the given input graph.
		
		Keyword arguments:
		G -- the input graph
		parameter -- a parameter value that determines the degree of sparsification
		attribute -- (optional) a previously calculated edge attribute. If none is provided, we will try to calculate it. 
		"""
		
		if attribute is None:
			attribute = self.getAttribute(G)
		
		return self._getBackbone(G, parameter, attribute)
	
	def getBackboneOfSize(self, G, edgeRatio, attribute=None):
		"""This is a convenience function that applies an appropriate parameterization
		algorithm to obtain a parameter value that yields a backbone of the desired size
		and then calls getBackbone(...) with that parameter value.
		
		Keyword arguments:
		G -- the input graph
		edgeRatio -- the target edge ratio
		attribute -- (optional) a previously calculated edge attribute. If none is provided, we will try to calculate it.
		"""
		if attribute is None:
			attribute = self.getAttribute(G)
		
		parameter = self.getParameter(G, edgeRatio, attribute)
		
		return self.getBackbone(G, parameter, attribute)
 
	@_abc.abstractmethod
	def _getBackbone(self, G, parameter, attribute):
		pass

	def getParameter(self, G, edgeRatio, attribute=None):
		""" This is a convenience function that applies an appropriate parameterization 
		algorithm to obtain a parameter value that yields a backbone of the desired size. """
		if attribute is None:
			attribute = self.getAttribute(G)

		paramAlgorithm = self._getParameterizationAlgorithm()
		parameter = paramAlgorithm.parameterize(self, G, attribute, edgeRatio)
		return parameter

	def _getParameterizationAlgorithm(self):
		return SimpleParameterization()

""" Represents an algorithm that, given a graph and a backbone algorithm, 
calculates a parameter value such that a desired edge ratio is met.""" 
class SimpleParameterization:
	def parameterize(self, algorithm, G, attribute, edgeRatio):
		return edgeRatio

""" Parameterizes a backbone algorithm using binary search """
class BinarySearchParameterization:
	
	def __init__(self, sizeIncreasesWithParameter, lowerParameterBound, upperParameterBound, maxSteps):
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
			backbone = algorithm._getBackbone(G, estimation, attribute)
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

""" Parameterizes a backbone algorithm using complete search (applicable only to algorithms which take as
input a parameter from a small set of possible values """
class CompleteSearchParameterization:
	
	def __init__(self, lowerParameterBound, upperParameterBound):
		self.lowerParameterBound = lowerParameterBound
		self.upperParameterBound = upperParameterBound
	
	def parameterize(self, algorithm, G, attribute, edgeRatio):
		bestParameter = self.lowerParameterBound
		bestRatio = 0.0
		minDistance = 100.0

		for i in range(self.lowerParameterBound, self.upperParameterBound + 1):
			backbone = algorithm._getBackbone(G, i, attribute)
			currentEdgeRatio = backbone.numberOfEdges() / G.numberOfEdges()

			distance = abs(currentEdgeRatio - edgeRatio)
			if distance < minDistance and abs(currentEdgeRatio) > _ABS_ZERO:
				minDistance = distance
				bestParameter = i
				bestRatio = currentEdgeRatio

		return bestParameter

""" Parametric variant of the Simmelian Backbones introduced by Nick et al. """
class SimmelianBackboneParametric(BackboneAlgorithm):

	def getAttribute(self, G):
		chiba = ChibaNishizekiTriangleCounter()
		triangles = chiba.getAttribute(G)
		so = SimmelianOverlapAttributizer(10)
		a_so = so.getAttribute(G, triangles)
		return a_so

	def _getBackbone(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, True)
		return gf.calculate(G, attribute)

	def _getParameterizationAlgorithm(self):
		return CompleteSearchParameterization(0, 10)
	
""" Non-parametric variant of the Simmelian Backbones introduced by Nick et al. """
class SimmelianBackboneNonParametric(BackboneAlgorithm):
	
	def getAttribute(self, G):
		chiba = ChibaNishizekiTriangleCounter()
		triangles = chiba.getAttribute(G)
		sj = SimmelianJaccardAttributizer()
		a_sj = sj.getAttribute(G, triangles)
		return a_sj
	
	def _getBackbone(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, True)
		return gf.calculate(G, attribute)

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

""" Multiscale Backbone that uses triangle counts as input edge weight. """
class SimmelianMultiscaleBackbone(BackboneAlgorithm):
	
	def getAttribute(self, G):
		chiba = ChibaNishizekiTriangleCounter()
		triangles = chiba.getAttribute(G)
		ms = MultiscaleAttributizer()
		a_ms = ms.getAttribute(G, triangles)
		return a_ms
       
	def _getBackbone(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)
	
	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)

""" An implementation of the Local Similarity sparsification approach introduced by Satuluri et al. """
class LocalSimilarityBackbone(BackboneAlgorithm):
	def getAttribute(self, G):
		attributizer = LocalSimilarityAttributizer()
		a_ls = attributizer.getAttribute(G, [])
		return a_ls

	def _getBackbone(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)

""" An implementation of the Multiscale backbone approach introduced by Serrano et al. """
class MultiscaleBackbone(BackboneAlgorithm):
	def getAttribute(self, G):
		inputAttribute = [0.0] * G.upperEdgeIdBound()
		for edge in G.edges():
			edgeId = G.edgeId(edge[0], edge[1])
			inputAttribute[edgeId] = G.weight(edge[0], edge[1])
		
		attributizer = MultiscaleAttributizer()
		attribute = attributizer.getAttribute(G, inputAttribute)
		return attribute

	def _getBackbone(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)
	
	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)
	

""" Random Edge sampling. Edges to keep in the backbone are selected uniformly at random. """
class RandomBackbone(BackboneAlgorithm):
	def getAttribute(self, G):
		attributizer = RandomAttributizer(1.0)
		a_r = attributizer.getAttribute(G, [1.0] * G.upperEdgeIdBound())
		return a_r
		
	def _getBackbone(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)
	
	def _getParameterizationAlgorithm(self):
		return SimpleParameterization()

""" A variant of the Forest Fire sparsification approach proposed by Leskovec et al. """
class ForestFireBackbone(BackboneAlgorithm):
	
	def __init__(self, burnProbability, targetBurntRatio):
		""" Creates a new instance of the Edge Forest Fire backbone algorithm.
		
		Keyword arguments:
		burnProbability -- the probability that the neighbor of a burnt node gets burnt as well.
		edgeRatio -- the fires will stop when a edgeRatio * edgeCount edges have been burnt.
		"""
		self.burnProbability = burnProbability
		self.targetBurntRatio = targetBurntRatio
	
	def getAttribute(self, G):
		attributizer = ForestFireAttributizer(self.burnProbability, self.targetBurntRatio)
		return attributizer.getAttribute(G, [])
		
	def _getBackbone(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, True)
		return gf.calculate(G, attribute)
	
	def _getParameterizationAlgorithm(self):
		 return BinarySearchParameterization(False, 0.0, 1.0, 20)

""" An implementation of the Local Degree backbone algorithm. """
class LocalDegreeBackbone(BackboneAlgorithm):
	
	def getAttribute(self, G):
		attributizer_ld = LocalDegreeAttributizer()
		a_ld = attributizer_ld.getAttribute(G, [])
		return a_ld
       
	def _getBackbone(self, G, parameter, attribute):
		gf = GlobalThresholdFilter(parameter, False)
		return gf.calculate(G, attribute)
	
	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)

	