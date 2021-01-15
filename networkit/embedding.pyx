# distutils: language=c++

from libc.stdint cimport uint64_t

from libcpp.vector cimport vector
from libcpp.string cimport string

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph

ctypedef uint64_t count

cdef extern from "<networkit/embedding/Node2Vec.hpp>":

    cdef cppclass _Node2Vec "NetworKit::Node2Vec"(_Algorithm):
        _Node2Vec(_Graph G, double P, double Q, count L, count N, count D) except +
        vector[vector[float]] getFeatures() except +

cdef class Node2Vec(Algorithm):
    """ 
    Algorithm to extract features from the graph with the node2vec(word2vec)
    algorithm according to [https://arxiv.org/pdf/1607.00653v1.pdf].
    CAUTION: 
    This algorithm could take a lot of time on large networks (many nodes).

    Node2Vec(G, P, Q, L, N, D)

    Create a Node2Vec algorithm object for Graph `G` with these params 
 
    Parameters
    ----------
    G : networkit.Graph
        The graph.
    P : double
        The ratio for returning to the previous node on a walk.
        P > max(Q,1) : less likely to sample an already-visited node in
                       the following two steps
        P < min(Q,1) : more likely to sample an already-visited node in
                       the following two steps
    Q : double
        The ratio for the direction of the next step
        Q > 1 : the random walk is biased towards nodes close to the
                previous one.
        Q < 1 : the random walk is biased towards nodes which are 
                further away from the previous one. 
    L : count
        The walk length.
    N : count
        The number of walks per node.
    D: count
        The dimension of the calculated embedding. 
    """

    cdef Graph _G
 
    def __cinit__(self, Graph G, P=1, Q=1, L=80, N=10, D=128):
        self._G = G
        self._this = new _Node2Vec(G._this, P, Q, L, N, D)

    def getFeatures(self):
        """
        Returns all feature vectors

        Returns
        -------
        A vector of vectors of floats.
        A vector containing feature vectors of all nodes
        """
        return (<_Node2Vec*>(self._this)).getFeatures()

