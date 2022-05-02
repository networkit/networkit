# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.map cimport map
from libcpp.vector cimport vector

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .matching cimport _Matching, Matching
from .structures cimport _Cover, Cover, _Partition, Partition, index, node

cdef extern from "<networkit/coarsening/GraphCoarsening.hpp>":

	cdef cppclass _GraphCoarsening "NetworKit::GraphCoarsening"(_Algorithm):
		_GraphCoarsening(_Graph) except +
		_Graph getCoarseGraph() except +
		vector[node] getFineToCoarseNodeMapping() except +
		map[node, vector[node]] getCoarseToFineNodeMapping() except +

cdef class GraphCoarsening(Algorithm):
	""" Abstract base class for coarsening measures"""

	cdef Graph _G

	def __init__(self, *args, **namedargs):
		if type(self) == GraphCoarsening:
			raise RuntimeError("Error, you may not use GraphCoarsening directly, use a sub-class instead")

	def getCoarseGraph(self):
		"""
		getCoarseGraph()

		Returns the coarse graph

		Returns
		-------
		networkit.Graph
			Coarse graph.
		"""
		return Graph(0).setThis((<_GraphCoarsening*>(self._this)).getCoarseGraph())

	def getFineToCoarseNodeMapping(self):
		"""
		getCoarseGraph()

		Returns the coarse graph

		Returns
		-------
		list(dict(int ``:`` int))
			List containing mapping, whereas index represents node u from fine graph and value represents node v from coarse graph  
		"""
		return (<_GraphCoarsening*>(self._this)).getFineToCoarseNodeMapping()

	def getCoarseToFineNodeMapping(self):
		"""
		getCoarseGraph()

		Returns the coarse graph

		Returns
		-------
		dict(int ``:`` list(int))
			Map containing node from coarse graph and all its mappings from finer graph
		"""
		return (<_GraphCoarsening*>(self._this)).getCoarseToFineNodeMapping()


cdef extern from "<networkit/coarsening/ParallelPartitionCoarsening.hpp>":

	cdef cppclass _ParallelPartitionCoarsening "NetworKit::ParallelPartitionCoarsening"(_GraphCoarsening):
		_ParallelPartitionCoarsening(_Graph, _Partition, bool_t) except +


cdef class ParallelPartitionCoarsening(GraphCoarsening):
	"""
	ParallelPartitionCoarsening(G, zeta, parallel = True)
	
	Coarsens graph according to a partition.
 	
 	Parameters
 	----------
 	G : networkit.Graph
		The input graph.
	zeta : networkit.Partition
		The partition, which is used for coarsening.
 	parallel : bool, optional
		If true, algorithm runs in parallel. Default: True
	"""
	def __cinit__(self, Graph G not None, Partition zeta not None, parallel = True):
		self._this = new _ParallelPartitionCoarsening(G._this, zeta._this, parallel)

cdef extern from "<networkit/coarsening/MatchingCoarsening.hpp>":

	cdef cppclass _MatchingCoarsening "NetworKit::MatchingCoarsening"(_GraphCoarsening):
		_MatchingCoarsening(_Graph, _Matching, bool_t) except +


cdef class MatchingCoarsening(GraphCoarsening):
	"""
	MatchingCoarsening(G, M, noSelfLoops=False)
	
	Coarsens graph according to a matching.
 	
 	Parameters
 	----------
 	G : networkit.Graph
		The input graph.
	M : networkit.matching.Matching
		The matching, which is used for coarsening.
 	noSelfLoops : bool, optional
		If true, self-loops are not produced. Default: False
	"""

	def __cinit__(self, Graph G not None, Matching M not None, bool_t noSelfLoops=False):
		self._this = new _MatchingCoarsening(G._this, M._this, noSelfLoops)
