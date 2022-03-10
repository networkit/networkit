# distutils: language=c++

from libc.stdint cimport uint64_t
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.map cimport map

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node

from .base cimport _Algorithm, Algorithm
from .dynamics cimport _GraphEvent, GraphEvent
from .graph cimport _Graph, Graph
from .structures cimport _Partition, Partition

cdef extern from "<networkit/components/ConnectedComponents.hpp>":

	cdef cppclass _ConnectedComponents "NetworKit::ConnectedComponents"(_Algorithm):
		_ConnectedComponents(_Graph G) except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		_Partition getPartition() except +
		map[index, count] getComponentSizes() except +
		vector[vector[node]] getComponents() except +
		@staticmethod
		_Graph extractLargestConnectedComponent(_Graph G, bool_t) nogil except +

cdef class ConnectedComponents(Algorithm):
	""" 
	ConnectedComponents(G)
	
	Determines the connected components and associated values for an undirected graph. 
	Create ConnectedComponents for Graph `G`.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G):
		self._G = G
		self._this = new _ConnectedComponents(G._this)

	def getPartition(self):
		""" 
		getPartition()

		Get a Partition that represents the components.

		Returns
		-------
		networkit.Partition
			A partition representing the found components.
		"""
		return Partition().setThis((<_ConnectedComponents*>(self._this)).getPartition())

	def numberOfComponents(self):
		""" 
		numberOfComponents()
		
		Get the number of connected components.

		Returns
		-------
		int
			The number of connected components.
		"""
		return (<_ConnectedComponents*>(self._this)).numberOfComponents()

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
		return (<_ConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponentSizes(self):
		""" 
		getComponentSizes()
		
		Get the component sizes.

		Returns
		-------
		dict(int ``:`` int)
			A dict containing the component ids and their size.
		"""
		return (<_ConnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		""" 
		getComponents()
		
		Get the connected components, each as a list of nodes.

		Returns
		-------
		list(int)
			The connected components.
		"""
		return (<_ConnectedComponents*>(self._this)).getComponents()

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
		graph : networkit.Graph
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

	cdef cppclass _ParallelConnectedComponents "NetworKit::ParallelConnectedComponents"(_Algorithm):
		_ParallelConnectedComponents(_Graph G, bool_t coarsening) except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		_Partition getPartition() except +
		vector[vector[node]] getComponents() except +


cdef class ParallelConnectedComponents(Algorithm):
	""" 
	ParallelConnectedComponents(G, coarsening=True)

	Determines the connected components and associated values for
	an undirected graph.

	Parameters
	----------
	graph : networkit.Graph
		The input graph
	coarsening : bool, optional
		Specifies whether the main algorithm based on label propagation (LP) 
		shall work recursively (true) or not (false) by coarsening/contracting 
		an LP-computed clustering. Defaults to true since we saw positive effects 
		in terms of running time for many networks. Beware of possible memory implications.
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G, coarsening=True	):
		self._G = G
		self._this = new _ParallelConnectedComponents(G._this, coarsening)

	def getPartition(self):
		""" 
		getPartition()
		
		Get a Partition that represents the components.

		Returns
		-------
		networkit.Partition
			A partition representing the found components.
		"""
		return Partition().setThis((<_ParallelConnectedComponents*>(self._this)).getPartition())

	def numberOfComponents(self):
		""" 
		numberOfComponents()
		
		Get the number of connected components.

		Returns
		-------
		int
			The number of connected components.
		"""
		return (<_ParallelConnectedComponents*>(self._this)).numberOfComponents()

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
			Component in which node `v` is situated
		"""
		return (<_ParallelConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponents(self):
		""" 
		getComponents()

		Get the connected components, each as a list of nodes.

		Returns
		-------
		list(int)
			The connected components.
		"""
		return (<_ParallelConnectedComponents*>(self._this)).getComponents()


cdef extern from "<networkit/components/StronglyConnectedComponents.hpp>":

	cdef cppclass _StronglyConnectedComponents "NetworKit::StronglyConnectedComponents":
		_StronglyConnectedComponents(_Graph G) except +
		void run() nogil except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		_Partition getPartition() except +
		map[node, count] getComponentSizes() except +
		vector[vector[node]] getComponents() except +


cdef class StronglyConnectedComponents:
	""" 
	StronglyConnectedComponents(G)

	Computes the strongly connected components of a directed graph.

	Parameters:
	-----------
	G : networkit.Graph
		The input graph.
	"""
	cdef _StronglyConnectedComponents* _this
	cdef Graph _G

	def __cinit__(self, Graph G):

		self._G = G
		self._this = new _StronglyConnectedComponents(G._this)

	def __dealloc__(self):
		del self._this

	def run(self):
		"""
		run()

		Runs the algorithm.
		"""
		with nogil:
			self._this.run()
		return self

	def getPartition(self):
		"""
		getPartition()

		Returns a Partition object representing the strongly connected components.

		Returns
		-------
		networkit.Partition
			The strongly connected components.
		"""
		return Partition().setThis(self._this.getPartition())

	def numberOfComponents(self):
		"""
		numberOfComponents()
		
		Returns the number of strongly connected components of the graph.

		Returns
		-------
		int
			The number of strongly connected components.
		"""
		return self._this.numberOfComponents()

	def componentOfNode(self, u):
		"""
		Returns the component of node `u`.

		Parameters
		----------
		u : int
			A node in the graph.

		Returns
		-------
		int
			The component of node `u`.
		"""
		return self._this.componentOfNode(u)

	def getComponentSizes(self):
		"""
		getComponentSizes()
		
		Returns a map with the component indexes as keys and their size as values.

		Returns
		-------
		dict(int ``:`` int)
			A dict with component indexes as keys and their size as values.
		"""
		return self._this.getComponentSizes()

	def getComponents(self):
		"""
		getComponents()
		
		Returns a list of components.

		Returns
		-------
		list(list(int))
			A list of components.
		"""
		return self._this.getComponents()


cdef extern from "<networkit/components/WeaklyConnectedComponents.hpp>":

	cdef cppclass _WeaklyConnectedComponents "NetworKit::WeaklyConnectedComponents"(_Algorithm):
		_WeaklyConnectedComponents(_Graph G) except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		map[index, count] getComponentSizes() except +
		vector[vector[node]] getComponents() except +

cdef class WeaklyConnectedComponents(Algorithm):
	""" 
	WeaklyConnectedComponents(G)

	Determines the weakly connected components of a directed graph.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _WeaklyConnectedComponents(G._this)

	def numberOfComponents(self):
		""" 
		numberOfComponents()

		Returns the number of components.

		Returns
		-------
		int
			The number of components.
		"""
		return (<_WeaklyConnectedComponents*>(self._this)).numberOfComponents()

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
			Component in which node `v` is situated
		"""
		return (<_WeaklyConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponentSizes(self):
		""" 
		getComponentSizes()
		
		Returns the map from component id to size.

		Returns
		-------
		dict(int ``:`` int)
			A dict that maps each component id to its size.
		"""
		return (<_WeaklyConnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		""" 
		getComponents()
		
		Returns all the components, each stored as (unordered) set of nodes.

		Returns
		-------
		list(list(int))
			A list of lists. Each inner vector contains all the nodes inside the component.
		"""
		return (<_WeaklyConnectedComponents*>(self._this)).getComponents()


cdef extern from "<networkit/components/BiconnectedComponents.hpp>":

	cdef cppclass _BiconnectedComponents "NetworKit::BiconnectedComponents"(_Algorithm):
		_BiconnectedComponents(_Graph G) except +
		count numberOfComponents() except +
		map[count, count] getComponentSizes() except +
		vector[vector[node]] getComponents() except +

cdef class BiconnectedComponents(Algorithm):
	""" 
	BiconnectedComponents(G)
	
	Determines the biconnected components of an undirected graph as defined in
	Tarjan, Robert. Depth-First Search and Linear Graph Algorithms. SIAM J.
	Comput. Vol 1, No. 2, June 1972.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
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
		dict(int ``:`` int)
			A dict that maps each component id to its size.
		"""
		return (<_BiconnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		""" 
		getComponents()
		
		Returns all the components, each stored as (unordered) set of nodes.

		Returns
		-------
		list(list(int))
			A list of lists. Each inner vector contains all the nodes inside the component.
		"""
		return (<_BiconnectedComponents*>(self._this)).getComponents()


cdef extern from "<networkit/components/DynConnectedComponents.hpp>":

	cdef cppclass _DynConnectedComponents "NetworKit::DynConnectedComponents"(_Algorithm):
		_DynConnectedComponents(_Graph G) except +
		void update(_GraphEvent) except +
		void updateBatch(vector[_GraphEvent]) except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		map[index, count] getComponentSizes() except +
		vector[vector[node]] getComponents() except +

cdef class DynConnectedComponents(Algorithm):
	""" 
	DynConnectedComponents(G)
	
	Determines and updates the connected components of an undirected graph.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _DynConnectedComponents(G._this)

	def numberOfComponents(self):
		"""
		numberOfComponents()
		
		Returns the number of components.

		Returns
		-------
		int
			The number of components.
		"""
		return (<_DynConnectedComponents*>(self._this)).numberOfComponents()

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
		return (<_DynConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponentSizes(self):
		""" 
		getComponentSizes()
		
		Returns the map from component id to size.

		Returns
		-------
		dict(int ``:`` int)
			A dict that maps each component id to its size.
		"""
		return (<_DynConnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		"""
		getComponents()
		
		Returns all the components, each stored as (unordered) set of nodes.

		Returns
		-------
		list(list(int))
			A list of lists. Each inner list contains all the nodes inside the component.
		"""
		return (<_DynConnectedComponents*>(self._this)).getComponents()

	def update(self, event):
		""" 
		update(event)
		
		Updates the connected components after an edge insertion or
		deletion.

		Parameters
		----------
		event : networkit.dynamics.GraphEvent
			The event that happened (edge deletion or insertion).
		"""
		(<_DynConnectedComponents*>(self._this)).update(_GraphEvent(event.type, event.u, event.v, event.w))

	def updateBatch(self, batch):
		""" 
		updateBatch(batch)
		
		Updates the connected components after a batch of edge insertions or
		deletions.

		Parameters
		----------
		batch : list(networkit.dynamics.GraphEvent)
			A list that contains a batch of edge insertions or deletions.
		"""
		cdef vector[_GraphEvent] _batch
		for event in batch:
			_batch.push_back(_GraphEvent(event.type, event.u, event.v, event.w))
		(<_DynConnectedComponents*>(self._this)).updateBatch(_batch)



cdef extern from "<networkit/components/DynWeaklyConnectedComponents.hpp>":

	cdef cppclass _DynWeaklyConnectedComponents "NetworKit::DynWeaklyConnectedComponents"(_Algorithm):
		_DynWeaklyConnectedComponents(_Graph G) except +
		void update(_GraphEvent) except +
		void updateBatch(vector[_GraphEvent]) except +
		count numberOfComponents() except +
		count componentOfNode(node query) except +
		map[index, count] getComponentSizes() except +
		vector[vector[node]] getComponents() except +

cdef class DynWeaklyConnectedComponents(Algorithm):
	""" 
	DynWeaklyConnectedComponents(G)
	
	Determines and updates the weakly connected components of a directed graph.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _DynWeaklyConnectedComponents(G._this)

	def numberOfComponents(self):
		""" 
		numberOfComponents()
		
		Returns the number of components.

		Returns
		-------
		int
			The number of components.
		"""
		return (<_DynWeaklyConnectedComponents*>(self._this)).numberOfComponents()

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
			Component in which node `v` is situated
		"""
		return (<_DynWeaklyConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponentSizes(self):
		""" 
		getComponentSizes()
		
		Returns the map from component id to size.

		Returns
		-------
		dict(int ``:`` int)
			A dict that maps each component id to its size.
		"""
		return (<_DynWeaklyConnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		""" 
		getComponents()
		
		Returns all the components, each stored as (unordered) set of nodes.

		Returns
		-------
		list(list(int))
			A list of lists. Each inner vector contains all the nodes
			inside the component.
		"""
		return (<_DynWeaklyConnectedComponents*>(self._this)).getComponents()

	def update(self, event):
		""" 
		update(event)
		
		Updates the connected components after an edge insertion or
		deletion.

		Parameters
		----------
		event : networkit.dynamics.GraphEvent
			The event that happened (edge deletion or insertion).
		"""
		(<_DynWeaklyConnectedComponents*>(self._this)).update(_GraphEvent(event.type, event.u, event.v, event.w))

	def updateBatch(self, batch):
		"""
		updateBatch(batch)
		
		Updates the connected components after a batch of edge insertions or
		deletions.

		Parameters
		----------
		batch : list(networkit.dynamics.GraphEvent)
			A vector that contains a batch of edge insertions or deletions.
		"""
		cdef vector[_GraphEvent] _batch
		for event in batch:
			_batch.push_back(_GraphEvent(event.type, event.u, event.v, event.w))
		(<_DynWeaklyConnectedComponents*>(self._this)).updateBatch(_batch)


