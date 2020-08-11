# distutils: language=c++

from libc.stdint cimport uint64_t
from libcpp cimport bool as bool_t

ctypedef uint64_t count

from .graph cimport _Graph, Graph


cdef extern from "<networkit/global/ClusteringCoefficient.hpp>" namespace "NetworKit::ClusteringCoefficient":

		double avgLocal(_Graph G, bool_t turbo) nogil except +
		double sequentialAvgLocal(_Graph G) nogil except +
		double approxAvgLocal(_Graph G, count trials) nogil except +
		double exactGlobal(_Graph G) nogil except +
		double approxGlobal(_Graph G, count trials) nogil except +

cdef class ClusteringCoefficient:
	@staticmethod
	def avgLocal(Graph G, bool_t turbo = False):
		"""
		DEPRECATED: Use centrality.LocalClusteringCoefficient and take average.

		This calculates the average local clustering coefficient of graph `G`. The graph may not contain self-loops.

		Parameters:
		-----------
		G : networkit.Graph
			The graph.

		Notes:
		------

		.. math:: c(G) := \\frac{1}{n} \sum_{u \in V} c(u)

		where

		.. math:: c(u) := \\frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}

		"""
		cdef double ret
		with nogil:
			ret = avgLocal(G._this, turbo)
		return ret

	@staticmethod
	def sequentialAvgLocal(Graph G):
		""" This calculates the average local clustering coefficient of graph `G` using inherently sequential triangle counting.
		Parameters:
		-----------
		G : networkit.Graph
			The graph.

		Notes:
		------

		.. math:: c(G) := \\frac{1}{n} \sum_{u \in V} c(u)

		where

		.. math:: c(u) := \\frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}

		"""
		cdef double ret
		with nogil:
			ret = sequentialAvgLocal(G._this)
		return ret

	@staticmethod
	def approxAvgLocal(Graph G, count trials):
		cdef double ret
		with nogil:
			ret = approxAvgLocal(G._this, trials)
		return ret

	@staticmethod
	def exactGlobal(Graph G):
		""" This calculates the global clustering coefficient. """
		cdef double ret
		with nogil:
			ret = exactGlobal(G._this)
		return ret

	@staticmethod
	def approxGlobal(Graph G, count trials):
		cdef double ret
		with nogil:
			ret = approxGlobal(G._this, trials)
		return ret

#external imports
import math
import logging

def clustering(G, error=0.01):
	"""
		Returns approximate average local clustering coefficient
		The maximum error can be given as a parameter and determines
		the number of samples taken.

		for details see:
			Schank, Wagner: Approximating Clustering Coefficient and Transitivity
	"""
	if G.numberOfNodes() < 100:
		return ClusteringCoefficient().avgLocal(G)
	else:
		nSamples = math.ceil(math.log(10) / (error**2)) # fixed confidence of 90%
		logging.info("taking {0} samples".format(nSamples))
		return ClusteringCoefficient().approxAvgLocal(G, nSamples)
