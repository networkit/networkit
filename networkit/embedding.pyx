# distutils: language=c++

from libcpp.vector cimport vector
from libcpp.string cimport string

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport count

cdef extern from "<networkit/embedding/Node2Vec.hpp>":

	cdef cppclass _Node2Vec "NetworKit::Node2Vec"(_Algorithm):
		_Node2Vec(_Graph G, double P, double Q, count L, count N, count D) except +
		vector[vector[float]] &getFeatures() except +

cdef class Node2Vec(Algorithm):
	""" 
	Node2Vec(G, P, Q, L, N, D)

	Algorithm to extract features from the graph with the node2vec(word2vec)
	algorithm according to [https://arxiv.org/pdf/1607.00653v1.pdf].

	Node2Vec learns embeddings for nodes in a graph by optimizing a neighborhood preserving
	objective. In order to achieve this, biased random walks are initiated for every node and the
	result is of probabilistic nature. Several input parameters control the specific behavior of
	the random walks. Amongst others Node2Vec is able to produce embeddings for visualization
	(D=2 or D=3) and machine learning (D=128 [default]). Both directed and undirected graphs
	withouth isolated nodes are supported.

	Note
	---- 
	This algorithm could take a lot of time on large networks (many nodes).
 
	Parameters
	----------
	G : networkit.Graph
		The graph.
	P : float
		The ratio for returning to the previous node on a walk.
		For P > max(Q,1) it is less likely to sample an already-visited node in the following two steps.
		For P < min(Q,1) it is more likely to sample an already-visited node in the following two steps.
	Q : float
		The ratio for the direction of the next step
		For Q > 1 the random walk is biased towards nodes close to the previous one.
		For Q < 1 the random walk is biased towards nodes which are further away from the previous one. 
	L : int
		The walk length.
	N : int
		The number of walks per node.
	D: int
		The dimension of the calculated embedding. 
	"""

	cdef Graph _G
 
	def __cinit__(self, Graph G, P=1, Q=1, L=80, N=10, D=128):
		self._G = G
		self._this = new _Node2Vec(G._this, P, Q, L, N, D)

	def getFeatures(self):
		"""
		getFeatures()

		Returns all feature vectors

		Returns
		-------
		list(list(float))
			A vector containing feature vectors of all nodes
		"""
		return (<_Node2Vec*>(self._this)).getFeatures()
