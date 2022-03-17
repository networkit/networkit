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
	
	Parameters
	----------
	G : networkit.graph
		The input graph.
	data : list(float)
		Numerical node value array.
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

