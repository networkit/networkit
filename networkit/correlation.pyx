# distutils: language=c++

from libcpp.vector cimport vector

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition

cdef extern from "<networkit/correlation/Assortativity.hpp>":

	cdef cppclass _Assortativity "NetworKit::Assortativity"(_Algorithm):
		_Assortativity(_Graph, vector[double]) except +
		_Assortativity(_Graph, _Partition) except +
		double getCoefficient() except +

cdef class Assortativity(Algorithm):
	""" 
	Assortativity(G, data)
	
	Assortativity computes a coefficient that expresses the correlation of a
	node attribute among connected pairs of nodes.
	
	For a partition of the node set (a nominal attribute), the coefficient is
	computed from the mixing matrix e, where e_ij is the fraction of edge counts
	from nodes in subset i to nodes in subset j. In undirected graphs, each edge
	{u, v} with u != v is counted once in each direction and each self-loop is
	counted once, and the mixing matrix is normalized by the total of these
	counts (2m minus the number of self-loops). In directed graphs, each edge is
	counted once and the matrix is normalized by m.
	
	Parameters
	----------
	G : networkit.graph
		The input graph.
	data : list(float) or networkit.Partition
		Numerical node value array, or a partition of the node set.
	"""
	cdef Graph G
	cdef vector[double] attribute
	cdef Partition partition

	def __cinit__(self, Graph G, data):
		if isinstance(data, Partition):
			self._this = new _Assortativity(G._this, (<Partition>data)._this)
			self.partition = <Partition>data
		else:
			self.attribute = <vector[double]?>data
			self._this = new _Assortativity(G._this, self.attribute)
		self.G = G

	def getCoefficient(self):
		"""
		getCoefficient()

		Return the assortativity coefficient.

		Returns
		-------
		float
			The assortativity coefficient.
		"""
		return (<_Assortativity*>(self._this)).getCoefficient()

