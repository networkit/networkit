""" This module contains algorithms for the sparsification of networks """

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

	@abstractmethod
	def getBackbone(self, graph, attribute, parameter):
		"""Returns a sparsified version of the given input graph.
		
		Keyword arguments:
		graph -- the input graph
		attribute -- a previously calculated edge attribute
		parameter -- a parameter value that determines the degree of sparsification 
		"""
		pass

	@abstractmethod 
	def getParameterizationAlgorithm(self):
		return SimpleParameterization()

	def parameterize(self, graph, attribute, edgeRatio):
		"""Returns a parameter value that will yield a backbone graph of approximately the desired size
		
		Keyword arguments:
		graph -- the input graph
		attribute -- a previously calculated edge attribute
		edgeRatio -- the target edge ratio
		"""
		pass

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

	def getBackbone(self, graph, attribute, parameter):
		gf = GlobalThresholdFilter(value, True)
		return gf.calculate(graph, attribute)

	def getParameterizationAlgorithm(self):
		return CompleteSearchParameterization(0, 10)
