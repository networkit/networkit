# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.unordered_set cimport unordered_set

from .base cimport _Algorithm, Algorithm
from .dynbase cimport _DynAlgorithm
from .dynbase import DynAlgorithm
from .dynamics cimport _GraphEvent, GraphEvent
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition, count, index, node

cdef extern from "<networkit/components/ComponentDecomposition.hpp>":
	cdef cppclass _ComponentDecomposition "NetworKit::ComponentDecomposition"(_Algorithm):
		_ComponentDecomposition(_Graph G) except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		_Partition getPartition() except +
		map[index, count] getComponentSizes() except +
		vector[vector[node]] getComponents() except +

cdef class ComponentDecomposition(Algorithm):
	cdef Graph _G

	def __init__(self, *args, **kwargs):
		if type(self) == ComponentDecomposition:
			raise RuntimeError("Error, you may not use ComponentDecomposition directly, use a sub-class instead")

	def getPartition(self):
		"""
		getPartition()

		Get a Partition that represents the components.

		Returns
		-------
		networkit.Partition
			A partition representing the found components.
		"""
		return Partition().setThis((<_ComponentDecomposition*>(self._this)).getPartition())

	def numberOfComponents(self):
		"""
		numberOfComponents()

		Get the number of connected components.

		Returns
		-------
		int
			The number of connected components.
		"""
		return (<_ComponentDecomposition*>(self._this)).numberOfComponents()

	def componentOfNode(self, v):
		"""
		componentOfNode(v)

		Get the the component in which node `v` is situated.

		Parameters
		----------
		v : int
			The node whose component is asked for.

		Returns
		-------
		int
			Component in which node `v` is situated.
		"""
		return (<_ComponentDecomposition*>(self._this)).componentOfNode(v)

	def getComponentSizes(self):
		"""
		getComponentSizes()

		Get the component sizes.

		Returns
		-------
		dict(int ``:`` int)
			A dict containing the component ids and their size.
		"""
		return (<_ComponentDecomposition*>(self._this)).getComponentSizes()

	def getComponents(self):
		"""
		getComponents()

		Get the connected components, each as a list of nodes.

		Returns
		-------
		list(int)
			The connected components.
		"""
		return (<_ComponentDecomposition*>(self._this)).getComponents()


cdef extern from "<networkit/components/ConnectedComponents.hpp>":

	cdef cppclass _ConnectedComponents "NetworKit::ConnectedComponents"(_ComponentDecomposition):
		_ConnectedComponents(_Graph G) except +
		@staticmethod
		_Graph extractLargestConnectedComponent(_Graph G, bool_t) except + nogil

cdef class ConnectedComponents(ComponentDecomposition):
	"""
	ConnectedComponents(G)

	Determines the connected components and associated values for an undirected graph.
	Create ConnectedComponents for Graph `G`.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	"""

	def __cinit__(self,  Graph G):
		self._this = new _ConnectedComponents(G._this)

	@staticmethod
	def extractLargestConnectedComponent(Graph graph, bool_t compactGraph = False):
		"""
		extractLargestConnectedComponent(graph, compactGraph=False)

		Constructs a new graph that contains only the nodes inside the
		largest connected component.

		Notes
		-----
		Available for undirected graphs only.

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		compactGraph : bool, optional
			If true, the node ids of the output graph will be compacted
			(i.e., re-numbered from 0 to n-1). If false, the node ids
			will not be changed. Default: False

		Returns
		-------
		networkit.Graph
			A graph that contains only the nodes inside the largest
			connected component.
		"""
		return Graph().setThis(_ConnectedComponents.extractLargestConnectedComponent(graph._this, compactGraph))

cdef extern from "<networkit/components/ParallelConnectedComponents.hpp>":

	cdef cppclass _ParallelConnectedComponents "NetworKit::ParallelConnectedComponents"(_ComponentDecomposition):
		_ParallelConnectedComponents(_Graph G, bool_t coarsening) except +

cdef class ParallelConnectedComponents(ComponentDecomposition):
	"""
	ParallelConnectedComponents(G, coarsening=True)

	Determines the connected components and associated values for
	an undirected graph.

	Parameters
	----------
	G : networkit.Graph
		The input graph
	coarsening : bool, optional
		Specifies whether the main algorithm based on label propagation (LP)
		shall work recursively (true) or not (false) by coarsening/contracting
		an LP-computed clustering. Defaults to true since we saw positive effects
		in terms of running time for many networks. Beware of possible memory implications. Default: True
	"""

	def __cinit__(self,  Graph G, coarsening=True	):
		self._this = new _ParallelConnectedComponents(G._this, coarsening)

cdef extern from "<networkit/components/StronglyConnectedComponents.hpp>":

	cdef cppclass _StronglyConnectedComponents "NetworKit::StronglyConnectedComponents"(_ComponentDecomposition):
		_StronglyConnectedComponents(_Graph G) except +

cdef class StronglyConnectedComponents(ComponentDecomposition):
	"""
	StronglyConnectedComponents(G)

	Computes the strongly connected components of a directed graph.

	Parameters:
	-----------
	G : networkit.Graph
		The input graph.
	"""

	def __cinit__(self, Graph G):
		self._this = new _StronglyConnectedComponents(G._this)

cdef extern from "<networkit/components/WeaklyConnectedComponents.hpp>":

	cdef cppclass _WeaklyConnectedComponents "NetworKit::WeaklyConnectedComponents"(_ComponentDecomposition):
		_WeaklyConnectedComponents(_Graph G) except +

cdef class WeaklyConnectedComponents(ComponentDecomposition):
	"""
	WeaklyConnectedComponents(G)

	Determines the weakly connected components of a directed graph.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	"""

	def __cinit__(self, Graph G):
		self._this = new _WeaklyConnectedComponents(G._this)

cdef extern from "<networkit/components/BiconnectedComponents.hpp>":

	cdef cppclass _BiconnectedComponents "NetworKit::BiconnectedComponents"(_Algorithm):
		_BiconnectedComponents(_Graph G) except +
		count numberOfComponents() except +
		map[index, count] getComponentSizes() except +
		vector[vector[count]] getComponents() except +
		unordered_set[node] getComponentsOfNode(node u) except +

cdef class BiconnectedComponents(Algorithm):
	"""
	BiconnectedComponents()

	Determines the biconnected components of an undirected graph as defined in
	Tarjan, Robert. Depth-First Search and Linear Graph Algorithms. SIAM J.
	Comput. Vol 1, No. 2, June 1972.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""

	def __cinit__(self, Graph G):
		self._this = new _BiconnectedComponents(G._this)

	def numberOfComponents(self):
		"""
		numberOfComponents()
		
		Returns the number of components.

		Returns
		-------
		int
			The number of components.
		"""
		return (<_BiconnectedComponents*>(self._this)).numberOfComponents()

	def getComponentSizes(self):
		"""
		getComponentSizes()
		
		Returns the map from component id to size.

		Returns
		-------
		dict(int : int)
			A dict that maps each component id to its size.
		"""
		return (<_BiconnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		"""
		getComponents()
		
		Returns all the components, each stored as (unordered) set of nodes.

		Returns
		-------
		vector[vector[node]]
			A vector of vectors. Each inner vector contains all the nodes inside the component.
		"""
		return (<_BiconnectedComponents*>(self._this)).getComponents()

	def getComponentsOfNode(self, node u):
		"""
		getComponentsOfNode(u)

		Get the components that contain node u.

		Parameters
		----------
		u : int
			The node.

		Returns
		-------
		set
			Components that contain node u.
		"""
		return (<_BiconnectedComponents*>(self._this)).getComponentsOfNode(u)

cdef extern from "<networkit/components/DynConnectedComponents.hpp>":

	cdef cppclass _DynConnectedComponents "NetworKit::DynConnectedComponents"(_ComponentDecomposition, _DynAlgorithm):
		_DynConnectedComponents(_Graph G) except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		map[index, count] getComponentSizes() except +
		vector[vector[node]] getComponents() except +

cdef class DynConnectedComponents(ComponentDecomposition, DynAlgorithm):
	"""
	DynConnectedComponents(G)

	Determines and updates the connected components of an undirected graph.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""

	def __cinit__(self, Graph G):
		self._this = new _DynConnectedComponents(G._this)



cdef extern from "<networkit/components/DynWeaklyConnectedComponents.hpp>":

	cdef cppclass _DynWeaklyConnectedComponents "NetworKit::DynWeaklyConnectedComponents"(_ComponentDecomposition, _DynAlgorithm):
		_DynWeaklyConnectedComponents(_Graph G) except +

cdef class DynWeaklyConnectedComponents(ComponentDecomposition, DynAlgorithm):
	"""
	DynWeaklyConnectedComponents(G)

	Determines and updates the weakly connected components of a directed graph.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""

	def __cinit__(self, Graph G):
		self._this = new _DynWeaklyConnectedComponents(G._this)

