# distutils: language=c++

from libcpp cimport bool as bool_t

from .graph cimport _Graph, Graph
from .structures cimport count


cdef extern from "<networkit/global/ClusteringCoefficient.hpp>" namespace "NetworKit::ClusteringCoefficient":

		double avgLocal(_Graph G, bool_t turbo) except + nogil
		double sequentialAvgLocal(_Graph G) except + nogil
		double approxAvgLocal(_Graph G, count trials) except + nogil
		double exactGlobal(_Graph G) except + nogil
		double approxGlobal(_Graph G, count trials) except + nogil

cdef class ClusteringCoefficient:
	"""
	Class, which provides static functions for computing additional information for clustering coefficients.
	A :code:`ClusteringCoefficient` object itself doesn't have to be created.
	"""
	@staticmethod
	def sequentialAvgLocal(Graph G):
		""" 
		sequentialAvgLocal(G)
		
		This calculates the average local clustering coefficient of graph `G` using inherently sequential triangle counting.

		Notes
		-----

		.. math:: c(G) := \\frac{1}{n} \sum_{u \in V} c(u)

		where

		.. math:: c(u) := \\frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		cdef double ret
		with nogil:
			ret = sequentialAvgLocal(G._this)
		return ret

	@staticmethod
	def approxAvgLocal(Graph G, count trials):
		"""
		approxAvgLocal(G, trials)

		Approximates the average local clustering coefficient.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		trials : int
			Number of runs. Higher values result in higher quality and larger running times.
		"""
		cdef double ret
		with nogil:
			ret = approxAvgLocal(G._this, trials)
		return ret

	@staticmethod
	def exactGlobal(Graph G):
		""" 
		exactGlobal(G)
		
		Calculates the global clustering coefficient.
		
		Parameters
		----------
		G : networkit.Graph
			The input graph.
		"""
		cdef double ret
		with nogil:
			ret = exactGlobal(G._this)
		return ret

	@staticmethod
	def approxGlobal(Graph G, count trials):
		"""
		approxGlobal(G, trials)

		Approximates the global clustering coefficient.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		trials : int
			Number of runs. Higher values result in higher quality and larger running times.
		"""
		cdef double ret
		with nogil:
			ret = approxGlobal(G._this, trials)
		return ret

#external imports
import math
import logging

def clustering(G, error=0.01):
	"""
	clustering(G, error=0.01)

	Returns approximate average local clustering coefficient. The maximum error can be given as a parameter 
	and determines the number of samples taken.

	For details see: Schank, Wagner: Approximating Clustering Coefficient and Transitivity

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	error : float
		Maximum allowed error. Default: 0.01
	"""
	if G.numberOfNodes() < 100:
		return ClusteringCoefficient().avgLocal(G)
	else:
		nSamples = math.ceil(math.log(10) / (error**2)) # fixed confidence of 90%
		logging.info("taking {0} samples".format(nSamples))
		return ClusteringCoefficient().approxAvgLocal(G, nSamples)
