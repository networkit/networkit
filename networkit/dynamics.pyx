# distutils: language=c++

from .base cimport _Algorithm
from .base cimport Algorithm
from .graph cimport _Graph, Graph
from .helpers import stdstring

def graphFromStream(stream, weighted, directed):
	""" 
	graphFromStream(stream, weighted, directed)

	Convenience function for creating a new graph from a stream of graph events

	Parameters
	----------
	stream : list(networkit.dynamics.GraphEvent)
		Event stream
	weighted : bool
		Produce a weighted or unweighted graph
	directed : bool
		Produce a directed or undirected graph
	"""
	G = Graph(0, weighted, directed)
	gu = GraphUpdater(G)
	gu.update(stream)
	return G

cdef class GraphEvent:
	"""
	GraphEvent(type, u, v, w)

	Representation of a graph event.

	Parameter :code:`type` is one of the following: 
	
	- networkit.dynamics.GraphEvent.NODE_ADDITION
	- networkit.dynamics.GraphEvent.NODE_REMOVAL
	- networkit.dynamics.GraphEvent.NODE_RESTORATION
	- networkit.dynamics.GraphEvent.EDGE_ADDITION
	- networkit.dynamics.GraphEvent.EDGE_REMOVAL
	- networkit.dynamics.GraphEvent.EDGE_WEIGHT_UPDATE
	- networkit.dynamics.GraphEvent.EDGE_WEIGHT_INCREMENT

	Parameters
	----------
	type: networkit.dynamics.GraphEvent.type
		Type of graph event.
	u : int
		Node u involved in graph event.
	v : int
		Node v involved in graph event.
	w : int, float
		Weight of edge between node u and v.
	"""
	NODE_ADDITION = 0
	NODE_REMOVAL = 1
	NODE_RESTORATION = 2
	EDGE_ADDITION = 3
	EDGE_REMOVAL = 4
	EDGE_WEIGHT_UPDATE = 5
	EDGE_WEIGHT_INCREMENT = 6
	TIME_STEP = 7

	@property
	def type(self):
		"""
		Property of networkit.dynamics.GraphEvent
		
		Type of graph event.
		"""
		return self._this.type
	
	@type.setter
	def type(self, t):
		self._this.type = t

	@property
	def u(self):
		"""
		Property of networkit.dynamics.GraphEvent
		
		Node u involved in graph event.
		"""
		return self._this.u
	
	@u.setter
	def u(self, u):
		self._this.u = u

	@property
	def v(self):
		"""
		Property of networkit.dynamics.GraphEvent
		
		Node v involved in graph event.
		"""
		return self._this.v
	
	@v.setter
	def v(self, v):
		self._this.v = v

	@property
	def w(self):
		"""
		Property of networkit.dynamics.GraphEvent
		
		Edgeweight w involved in graph event.
		"""
		return self._this.w
	
	@w.setter
	def w(self, w):
		self._this.w = w

	def __cinit__(self, _GraphEventType type, node u, node v, edgeweight w):
		self._this = _GraphEvent(type, u, v, w)

	def toString(self):
		return self._this.toString().decode("utf-8")

	def __repr__(self):
		return self.toString()

	def __eq__(self, GraphEvent other not None):
		return self._this == other._this

	def __ne__(self, GraphEvent other not None):
		return self._this != other._this

	def __lt__(self, GraphEvent other not None):
		return self._this < other._this

	def __gt__(self, GraphEvent other not None):
		return self._this > other._this

	def __le__(self, GraphEvent other not None):
		return self._this <= other._this

	def __ge__(self, GraphEvent other not None):
		return self._this >= other._this

cdef extern from "<networkit/dynamics/DGSStreamParser.hpp>":

	cdef cppclass _DGSStreamParser "NetworKit::DGSStreamParser":
		_DGSStreamParser(string path, bool_t mapped, node baseIndex) except +
		vector[_GraphEvent] getStream() except +

cdef class DGSStreamParser:
	"""
	DGSStreamParser(path, mapped=True, baseIndex=0)

	Create a DGSStreamParser, handling graph streams encoded in DGS-based files.
	For documentation about DGS, see: https://graphstream-project.org/doc/Advanced-Concepts/The-DGS-File-Format/

	Parameters
	----------
	path : str
		Filename including path for DGS-file.
	mapped : bool, optional
		Indicates whether file includes mapping. Default: True
	baseIndex : int, optional
		Indicates whether a fixed base index should be set. Default: 0
	"""
	cdef _DGSStreamParser* _this

	def __cinit__(self, path, mapped=True, baseIndex=0):
		self._this = new _DGSStreamParser(stdstring(path), mapped, baseIndex)

	def __dealloc__(self):
		del self._this

	def getStream(self):
		"""
		getStream()

		Returns a list of graph events (networkit.GraphEvent).

		Returns
		-------
		list(networkit.dynamics.GraphEvent)
			A list of graph events.
		"""

		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.getStream()]

cdef extern from "<networkit/dynamics/DGSWriter.hpp>":

	cdef cppclass _DGSWriter "NetworKit::DGSWriter":
		void write(vector[_GraphEvent] stream, string path) except +


cdef class DGSWriter:
	"""
	DGSWriter()

	Creates a DGS writer object.
	"""
	cdef _DGSWriter* _this

	def __cinit__(self):
		self._this = new _DGSWriter()

	def __dealloc__(self):
		del self._this

	def write(self, stream, path):
		"""
		write(stream, path)

		Stores a graph event stream in a file (DGS format).

		Parameters
		----------
		stream : list(networkit.dynamics.GraphEvent)
			A list of graph events.
		path : str
			File to write the stream to.
		"""
		cdef vector[_GraphEvent] _stream
		for ev in stream:
			_stream.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.write(_stream, stdstring(path))


cdef extern from "<networkit/dynamics/GraphDifference.hpp>":

	cdef cppclass _GraphDifference "NetworKit::GraphDifference"(_Algorithm):
		_GraphDifference(const _Graph &G1, const _Graph &G2) except +
		vector[_GraphEvent] &getEdits() except +
		count getNumberOfEdits() except +
		count getNumberOfNodeAdditions() except +
		count getNumberOfNodeRemovals() except +
		count getNumberOfNodeRestorations() except +
		count getNumberOfEdgeAdditions() except +
		count getNumberOfEdgeRemovals() except +
		count getNumberOfEdgeWeightUpdates() except +

cdef class GraphDifference(Algorithm):
	"""
	GraphDifference(G1, G2)

	Calculate the edge difference between two graphs.

	This calculates which graph edge additions or edge removals are
	necessary to transform one given graph into another given graph.

	Both graphs need to have the same node set, directed graphs are not
	supported currently.

	Note that edge weight differences are not detected but edge
	addition events set the correct edge weight.

	Parameters
	----------
	G1 : networkit.Graph
		The first graph to compare
	G2 : networkit.Graph
		The second graph to compare
	"""
	cdef Graph _G1, _G2

	def __cinit__(self, Graph G1, Graph G2):
		self._this = new _GraphDifference(G1._this, G2._this)
		self._G1 = G1
		self._G2 = G2

	def getEdits(self):
		""" 
		getEdits()

		Get the required edits.

		Returns
		-------
		list(networkit.dynamics.GraphEvent)
			A list of graph events.
		"""
		cdef _GraphEvent ev
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in (<_GraphDifference*>(self._this)).getEdits()]

	def getNumberOfEdits(self):
		""" 
		getNumberOfEdits()
		
		Get the required number of edits.

		Returns
		-------
		int
			The number of edits.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfEdits()

	def getNumberOfNodeAdditions(self):
		""" 
		getNumberOfNodeAdditions()

		Get the required number of node additions.

		Returns
		-------
		int
			The number of node additions.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfNodeAdditions()

	def getNumberOfNodeRemovals(self):
		""" 
		getNumberOfNodeRemovals()

		Get the required number of node removals.

		Returns
		-------
		int
			The number of node removals.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfNodeRemovals()

	def getNumberOfNodeRestorations(self):
		""" 
		getNumberOfNodeRestorations()

		Get the required number of node restorations.

		Returns
		-------
		int
			The number of node restorations.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfNodeRestorations()

	def getNumberOfEdgeAdditions(self):
		""" 
		getNumberOfEdgeAdditions()
		
		Get the required number of edge additions.

		Returns
		-------
		int
			The number of edge additions.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfEdgeAdditions()

	def getNumberOfEdgeRemovals(self):
		""" 
		getNumberOfEdgeRemovals()
		
		Get the required number of edge removals.

		Returns
		-------
		int
			The number of edge removals.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfEdgeRemovals()

	def getNumberOfEdgeWeightUpdates(self):
		""" 
		getNumberOfEdgeWeightUpdates()
		
		Get the required number of edge weight updates.

		Returns
		-------
		int
			The number of edge weight updates.
		"""
		return (<_GraphDifference*>(self._this)).getNumberOfEdgeWeightUpdates()


cdef extern from "<networkit/dynamics/GraphUpdater.hpp>":

	cdef cppclass _GraphUpdater "NetworKit::GraphUpdater":
		_GraphUpdater(_Graph G) except +
		void update(vector[_GraphEvent] stream) nogil except +
		vector[pair[count, count]] &getSizeTimeline() except +

cdef class GraphUpdater:
	""" 
	GraphUpdater(G)
	
	Updates a graph according to a stream of graph events.

	Parameters
	----------
	G : networkit.Graph
		 Initial graph
	"""
	cdef _GraphUpdater* _this
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _GraphUpdater(G._this)

	def __dealloc__(self):
		del self._this

	def update(self, stream):
		"""
		update(stream)

		Update the underlying graph based on an input stream.

		Parameters
		----------
		stream : list(networkit.dynamics.GraphEvent)
			A list of graph events.
		"""

		cdef vector[_GraphEvent] _stream
		for ev in stream:
			_stream.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		with nogil:
			self._this.update(_stream)
