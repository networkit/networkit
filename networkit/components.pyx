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
	""" Determines the connected components and associated values for an undirected graph.

	ConnectedComponents(G)

	Create ConnectedComponents for Graph `G`.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G):
		self._G = G
		self._this = new _ConnectedComponents(G._this)

	def getPartition(self):
		""" Get a Partition that represents the components.

		Returns:
		--------
		networkit.Partition
			A partition representing the found components.
		"""
		return Partition().setThis((<_ConnectedComponents*>(self._this)).getPartition())

	def numberOfComponents(self):
		""" Get the number of connected components.

		Returns:
		--------
		count:
			The number of connected components.
		"""
		return (<_ConnectedComponents*>(self._this)).numberOfComponents()

	def componentOfNode(self, v):
		"""  Get the the component in which node `v` is situated.

		v : node
			The node whose component is asked for.
		"""
		return (<_ConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponentSizes(self):
		""" Get the component sizes.

		Returns:
		--------
		map:
			The map from component to size.
		"""
		return (<_ConnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		""" Get the connected components, each as a list of nodes.

		Returns:
		--------
		list:
			The connected components.
		"""
		return (<_ConnectedComponents*>(self._this)).getComponents()

	@staticmethod
	def extractLargestConnectedComponent(Graph graph, bool_t compactGraph = False):
		"""
			Constructs a new graph that contains only the nodes inside the
			largest connected component.

			Parameters:
			-----------
			graph: networkit.Graph
				The input graph
			compactGraph: bool
				if true, the node ids of the output graph will be compacted
				(i.e., re-numbered from 0 to n-1). If false, the node ids
				will not be changed.

			Returns:
			--------
			networkit.Graph
				A graph that contains only the nodes inside the largest
				connected component.


			Note:
			-----
			Available for undirected graphs only.
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
	""" Determines the connected components and associated values for
		an undirected graph.
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G, coarsening=True	):
		self._G = G
		self._this = new _ParallelConnectedComponents(G._this, coarsening)

	def getPartition(self):
		""" Get a Partition that represents the components.

		Returns:
		--------
		networkit.Partition
			A partition representing the found components.
		"""
		return Partition().setThis((<_ParallelConnectedComponents*>(self._this)).getPartition())

	def numberOfComponents(self):
		""" Get the number of connected components.

		Returns:
		--------
		count:
			The number of connected components.
		"""
		return (<_ParallelConnectedComponents*>(self._this)).numberOfComponents()

	def componentOfNode(self, v):
		"""  Get the the component in which node `v` is situated.

		v : node
			The node whose component is asked for.
		"""
		return (<_ParallelConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponents(self):
		""" Get the connected components, each as a list of nodes.

		Returns:
		--------
		list:
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
	cdef _StronglyConnectedComponents* _this
	cdef Graph _G

	def __cinit__(self, Graph G):
		""" Computes the strongly connected components of a directed graph.

			Parameters:
			-----------
			G : networkit.Graph
				The graph.
		"""
		self._G = G
		self._this = new _StronglyConnectedComponents(G._this)

	def __dealloc__(self):
		del self._this

	def run(self):
		"""
		Runs the algorithm.
		"""
		with nogil:
			self._this.run()
		return self

	def getPartition(self):
		"""
		Returns a Partition object representing the strongly connected components.

		Returns:
		--------
		Partition
			The strongly connected components.
		"""
		return Partition().setThis(self._this.getPartition())

	def numberOfComponents(self):
		"""
		Returns the number of strongly connected components of the graph.

		Returns:
		--------
		int
			The number of strongly connected components.
		"""
		return self._this.numberOfComponents()

	def componentOfNode(self, u):
		"""
		Returns the component of node `u`.

		Parameters:
		-----------
		u : node
			A node in the graph.

		Returns:
		int
			The component of node `u`.
		"""
		return self._this.componentOfNode(u)

	def getComponentSizes(self):
		"""
		Returns a map with the component indexes as keys, and their size as values.

		Returns:
		--------
		map[index, count]
			Map with component indexes as keys, and their size as values.
		"""
		return self._this.getComponentSizes()

	def getComponents(self):
		"""
		Returns a list of components.

		Returns:
		--------
		list[list[node]]
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
	""" Determines the weakly connected components of a directed graph.

		Parameters:
		-----------
		G : networkit.Graph
			The graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _WeaklyConnectedComponents(G._this)

	def numberOfComponents(self):
		""" Returns the number of components.

			Returns:
			--------
			count
				The number of components.
		"""
		return (<_WeaklyConnectedComponents*>(self._this)).numberOfComponents()

	def componentOfNode(self, v):
		""" Returns the the component in which node @a u is.

			Parameters:
			-----------
			v : node
				The node.
		"""
		return (<_WeaklyConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponentSizes(self):
		""" Returns the map from component to size.

			Returns:
			--------
			map[index, count]
			 	A map that maps each component to its size.
		"""
		return (<_WeaklyConnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		""" Returns all the components, each stored as (unordered) set of nodes.

			Returns:
			--------
			vector[vector[node]]
				A vector of vectors. Each inner vector contains all the nodes inside the component.
		"""
		return (<_WeaklyConnectedComponents*>(self._this)).getComponents()


cdef extern from "<networkit/components/BiconnectedComponents.hpp>":

	cdef cppclass _BiconnectedComponents "NetworKit::BiconnectedComponents"(_Algorithm):
		_BiconnectedComponents(_Graph G) except +
		count numberOfComponents() except +
		map[count, count] getComponentSizes() except +
		vector[vector[node]] getComponents() except +

cdef class BiconnectedComponents(Algorithm):
	""" Determines the biconnected components of an undirected graph as defined in
		Tarjan, Robert. Depth-First Search and Linear Graph Algorithms. SIAM J.
		Comput. Vol 1, No. 2, June 1972.


		Parameters:
		-----------
		G : networkit.Graph
			The graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _BiconnectedComponents(G._this)

	def numberOfComponents(self):
		""" Returns the number of components.

			Returns:
			--------
			count
				The number of components.
		"""
		return (<_BiconnectedComponents*>(self._this)).numberOfComponents()

	def getComponentSizes(self):
		""" Returns the map from component to size.

			Returns:
			--------
			map[count, count]
			A map that maps each component to its size.
		"""
		return (<_BiconnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		""" Returns all the components, each stored as (unordered) set of nodes.

			Returns:
			--------
			vector[vector[node]]
				A vector of vectors. Each inner vector contains all the nodes inside the component.
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
	""" Determines and updates the connected components of an undirected graph.

		Parameters:
		-----------
		G : networkit.Graph
			The graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _DynConnectedComponents(G._this)

	def numberOfComponents(self):
		""" Returns the number of components.

			Returns:
			--------
			count
				The number of components.
		"""
		return (<_DynConnectedComponents*>(self._this)).numberOfComponents()

	def componentOfNode(self, v):
		""" Returns the the component in which node @a u is.

			Parameters:
			-----------
			v : node
				The node.
		"""
		return (<_DynConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponentSizes(self):
		""" Returns the map from component to size.

			Returns:
			--------
			map[index, count]
			 	A map that maps each component to its size.
		"""
		return (<_DynConnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		""" Returns all the components, each stored as (unordered) set of nodes.

			Returns:
			--------
			vector[vector[node]]
				A vector of vectors. Each inner vector contains all the nodes inside the component.
		"""
		return (<_DynConnectedComponents*>(self._this)).getComponents()

	def update(self, event):
		""" Updates the connected components after an edge insertion or
			deletion.

			Parameters:
			-----------
			event : GraphEvent
				The event that happened (edge deletion or insertion).
		"""
		(<_DynConnectedComponents*>(self._this)).update(_GraphEvent(event.type, event.u, event.v, event.w))

	def updateBatch(self, batch):
		""" Updates the connected components after a batch of edge insertions or
			deletions.

			Parameters:
			-----------
			batch : vector[GraphEvent]
				A vector that contains a batch of edge insertions or deletions.
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
	""" Determines and updates the weakly connected components of a directed graph.

		Parameters:
		-----------
		G : networkit.Graph
			The graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _DynWeaklyConnectedComponents(G._this)

	def numberOfComponents(self):
		""" Returns the number of components.

			Returns:
			--------
			count
				The number of components.
		"""
		return (<_DynWeaklyConnectedComponents*>(self._this)).numberOfComponents()

	def componentOfNode(self, v):
		""" Returns the the component in which node @a u is.

			Parameters:
			-----------
			v : node
				The node.
		"""
		return (<_DynWeaklyConnectedComponents*>(self._this)).componentOfNode(v)

	def getComponentSizes(self):
		""" Returns the map from component to size.

			Returns:
			--------
			map[index, count]
			 	A map that maps each component to its size.
		"""
		return (<_DynWeaklyConnectedComponents*>(self._this)).getComponentSizes()

	def getComponents(self):
		""" Returns all the components, each stored as (unordered) set of nodes.

			Returns:
			--------
			vector[vector[node]]
				A vector of vectors. Each inner vector contains all the nodes
				inside the component.

		"""
		return (<_DynWeaklyConnectedComponents*>(self._this)).getComponents()

	def update(self, event):
		""" Updates the connected components after an edge insertion or
			deletion.

			Parameters:
			-----------
			event : GraphEvent
				The event that happened (edge deletion or insertion).
		"""
		(<_DynWeaklyConnectedComponents*>(self._this)).update(_GraphEvent(event.type, event.u, event.v, event.w))

	def updateBatch(self, batch):
		""" Updates the connected components after a batch of edge insertions or
			deletions.

			Parameters:
			-----------
			batch : vector[GraphEvent]
				A vector that contains a batch of edge insertions or deletions.
		"""
		cdef vector[_GraphEvent] _batch
		for event in batch:
			_batch.push_back(_GraphEvent(event.type, event.u, event.v, event.w))
		(<_DynWeaklyConnectedComponents*>(self._this)).updateBatch(_batch)


