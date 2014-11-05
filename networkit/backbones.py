""" This module contains algorithms for the sparsification of networks, i.e. Backbone algorithms. """

__author__ = "Gerd Lindner"

from abc import ABCMeta, abstractmethod
from _NetworKit import ChibaNishizekiTriangleCounter, SimmelianJaccardAttributizer, GlobalThresholdFilter, LocalSimilarityAttributizer, MultiscaleAttributizer, SimmelianOverlapAttributizer, RandomAttributizer, LocalDegreeAttributizer, ForestFireAttributizer

#SimmelianBackboneParametric, SimmelianBackboneNonParametric, MultiscaleBackbone, LocalSimilarityBackbone, SimmelianMultiscaleBackbone, RandomBackbone

""" Abstract base class representing a graph sparsification algorithm. 
 Sparsification algorithms can be split up into a precalculation step
 (calculation of an attribute) and a global filtering step. """
class BackboneAlgorithm(object):
	__metaclass__ = ABCMeta
	
	@abstractmethod
	def getAttribute(self, graph):
		"""Returns an edge attribute"""
		pass

	def getBackbone(self, graph, parameter, attribute=None):
		"""Returns a sparsified version of the given input graph.
		
		Keyword arguments:
		graph -- the input graph
		parameter -- a parameter value that determines the degree of sparsification
		attribute -- (optional) a previously calculated edge attribute. If none is provided, we will try to calculate it. 
		"""
		
		if attribute is None:
			attribute = self.getAttribute()
		
		return self._getBackbone(graph)
	
	def getBackboneOfSize(self, graph, edgeRatio, attribute=None):
		"""This is a convenience function that applies an appropriate parameterization
		algorithm to obtain a parameter value that yields a backbone of the desired size
		and then calls getBackbone(...) with that parameter value.
		
		Keyword arguments:
		graph -- the input graph
		edgeRatio -- the target edge ratio
		attribute -- (optional) a previously calculated edge attribute. If none is provided, we will try to calculate it.
		"""
		if attribute is None:
			attribute = self.getAttribute(graph)
		
		paramAlgorithm = _getParameterizationAlgorithm()
		parameter = paramAlgorithm.parameterize(graph, attribute, edgeRatio)
		
		return self.getBackbone(graph, parameter, attribute)
 
	@abstractmethod
	def _getBackbone(self, graph, parameter, attribute):
		pass
 
	def _getParameterizationAlgorithm(self):
		return SimpleParameterization()

""" Represents an algorithm that, given a graph and a backbone algorithm, 
calculates a parameter value such that a desired edge ratio is met.""" 
class SimpleParameterization:
	def parameterize(self, graph, attribute, edgeRatio):
		return edgeRatio

""" Parameterizes a backbone algorithm using binary search """
class BinarySearchParameterization:
	
	ABS_ZERO = 1e-5
	
	def __init__(self, sizeIncreasesWithParameter, lowerParameterBound, upperParameterBound, maxSteps):
		self.sizeIncreasesWithParameter = sizeIncreasesWithParameter
		self.lowerParameterBound = lowerParameterBound
		self.upperParameterBound = upperParameterBound
		self.maxSteps = maxSteps
	
	def parameterize(self, algorithm, graph, attribute, edgeRatio):
		lowerBound = self.lowerParameterBound
		upperBound = self.upperParameterBound
		estimation = self.lowerParameterBound
		bestParameter = self.lowerParameterBound
		minDistance = self.upperParameterBound

		for i in range(0, maxSteps):
			estimation = (lowerBound + upperBound) / 2.0
			backbone = algorithm.getBackbone(graph, attribute, estimation)
			currentEdgeRatio = backbone.numberOfEdges() / graph.numberOfEdges()

			distance = abs(currentEdgeRatio - targetEdgeRatio)

			if distance < minDistance and abs(currentEdgeRatio) > ABS_ZERO:
				minDistance = distance
				bestParameter = estimation

				#"Exact" hit?
				if abs(currentEdgeRatio - targetEdgeRatio) < 1e-7:
					break;

			increase = (currentEdgeRatio < targetEdgeRatio)
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
	
	ABS_ZERO = 1e-5
	
	def __init__(self, lowerParameterBound, upperParameterBound):
		self.lowerParameterBound = lowerParameterBound
		self.upperParameterBound = upperParameterBound
	
	def parameterize(self, algorithm, graph, attribute, edgeRatio):
		bestParameter = self.lowerParameterBound
		bestRatio = 0.0
		minDistance = 100.0

		for i in range(self.lowerParameterBound, self.upperParameterBound + 1):
			backbone = algorithm.getBackbone(graph, attribute, i)
			currentEdgeRatio = backbone.numberOfEdges() / graph.numberOfEdges()

			distance = abs(currentEdgeRatio - targetEdgeRatio)
			if distance < minDistance and abs(currentEdgeRatio) > ABS_ZERO:
				minDistance = distance
				
			bestParameter = i
			bestRatio = currentEdgeRatio

			return bestParameter

""" Parametric variant of the Simmelian Backbones introduced by Nick et al. """
class SimmelianBackboneParametric:

	def getAttribute(self, graph):
		chiba = ChibaNishizekiTriangleCounter()
		triangles = chiba.getAttribute(graph)
		sj = SimmelianJaccardAttributizer()
		a_sj = sj.getAttribute(graph, triangles)
		return a_sj

	def _getBackbone(self, graph, attribute, parameter):
		gf = GlobalThresholdFilter(value, True)
		return gf.calculate(graph, attribute)

	def _getParameterizationAlgorithm(self):
		return CompleteSearchParameterization(0, 10)
	
""" Non-parametric variant of the Simmelian Backbones introduced by Nick et al. """
class SimmelianBackboneNonParametric:
	
	def getAttribute(self, graph):
		chiba = backbones.ChibaNishizekiTriangleCounter()
		triangles = chiba.getAttribute(graph)
		sj = backbones.SimmelianJaccardAttributizer()
		a_sj = sj.getAttribute(graph, triangles)
		return a_sj
	
	def _getBackbone(self, graph, attribute, parameter):
		gf = backbones.GlobalThresholdFilter(value, True)
		return gf.calculate(graph, attribute)

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(False, 0.0, 1.0, 20)

""" Multiscale Backbone that uses triangle counts as input edge weight. """
class SimmelianMultiscaleBackbone:
	
	def getAttribute(self, graph):
		chiba = backbones.ChibaNishizekiTriangleCounter()
		triangles = chiba.getAttribute(graph)
		ms = backbones.MultiscaleAttributizer()
		a_ms = ms.getAttribute(graph, triangles)
		return a_ms
       
	def _getBackbone(self, graph, attribute, parameter):
		gf = backbones.GlobalThresholdFilter(value, False)
		return gf.calculate(graph, attribute)
	
	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)

""" An implementation of the Local Similarity sparsification approach introduced by Satuluri et al. """
class LocalSimilarityBackbone:
	def getAttribute(self, graph):
		attributizer = backbones.LocalSimilarityAttributizer()
		a_ls = attributizer.getAttribute(graph, [])
		return a_ls

	def _getBackbone(self, graph, attribute, parameter):
		gf = backbones.GlobalThresholdFilter(value, False)
		return gf.calculate(graph, attribute)

	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)

""" An implementation of the Multiscale backbone approach introduced by Serrano et al. """
class MultiscaleBackbone:
	def getAttribute(self, graph):
		# TODO we might use a precalculated edge attribute for speedup, but that
        # requires writable edge attributes in python.
		return None

	def _getBackbone(self, graph, attribute, parameter):
		msb = backbones.MultiscaleBackbone(value)
		return msb.calculate(graph)
	
	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)
	

""" Random Edge sampling. Edges to keep in the backbone are selected uniformly at random. """
class RandomBackbone:
	def getAttribute(self, graph):
		attributizer = backbones.RandomAttributizer()
		a_r = attributizer.getAttribute(graph, [1.0] * graph.upperEdgeIdBound())
		
	def _getBackbone(self, graph, attribute, parameter):
		gf = backbones.GlobalThresholdFilter(value, False)
	
	def _getParameterizationAlgorithm(self):
		return SimpleParameterization()

""" A variant of the Forest Fire sparsification approach proposed by Leskovec et al. """
class ForestFireBackbone:
	
	def __init__(self, burnProbability, targetBurntRatio):
		""" Creates a new instance of the Edge Forest Fire backbone algorithm.
		
		Keyword arguments:
		burnProbability -- the probability that the neighbor of a burnt node gets burnt as well.
		edgeRatio -- the fires will stop when a edgeRatio * edgeCount edges have been burnt.
		"""
		self.burnProbability = burnProbability
		self.targetBurntRatio = targetBurntRatio
	
	def getAttribute(self, graph):
		attributizer = backbones.ForestFireAttributizer(self.burnProbability, self.targetBurntRatio)
		return attributizer.getAttribute(graph, [])
		
	def _getBackbone(self, graph, attribute, parameter):
		gf = backbones.GlobalThresholdFilter(parameter, True)
		return gf.calculate(graph, attribute)
	
	def _getParameterizationAlgorithm(self):
		 return BinarySearchParameterization(False, 0.0, 1.0, 20)

""" An implementation of the Local Degree backbone algorithm. """
class LocalDegreeBackbone:
	
	def getAttribute(self, graph):
		attributizer_ld = backbones.LocalDegreeAttributizer()
		a_ld = attributizer_ld.getAttribute(graph, [])
		return a_ld
       
	def _getBackbone(self, graph, attribute, parameter):
		gf = backbones.GlobalThresholdFilter(parameter, False)
		return gf.calculate(graph, attribute)
	
	def _getParameterizationAlgorithm(self):
		return BinarySearchParameterization(True, 0.0, 1.0, 20)

	