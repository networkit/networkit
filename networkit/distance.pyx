# distutils: language=c++

from cython.operator import dereference, preincrement

from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.map cimport map
from libcpp cimport bool as bool_t
from libcpp.string cimport string
from libcpp.set cimport set
from libcpp.unordered_map cimport unordered_map

from .base cimport _Algorithm, Algorithm
from .dynbase cimport _DynAlgorithm
from .dynbase import DynAlgorithm
from .dynamics cimport _GraphEvent
from .graph cimport _Graph, Graph
from .helpers import stdstring
from .helpers cimport maybe_asarray_1d, maybe_asarray_2d
from .structures cimport count, index, node, edgeweight

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

none = _none

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef extern from "<networkit/distance/STSP.hpp>":
	cdef cppclass _STSP "NetworKit::STSP"(_Algorithm):
		_STSP(_Graph G, node source, node target, bool_t storePred) except +
		_STSP(_Graph G, node source, vector[node].iterator targetsFirst,
			vector[node].iterator targetsLast, bool_t storePred) except +
		void setSource(node newSource) except +
		void setTarget(node newTarget) except +
		void setTargets(vector[node].iterator targetsFirst, vector[node].iterator targetsLast) except +
		vector[node] getPath() except +
		vector[node] getPredecessors() except +
		edgeweight getDistance() except +
		vector[edgeweight] &getDistances() except +
		unordered_map[node, index] &getTargetIndexMap() except +

cdef class STSP(Algorithm):
	""" 
	STSP(G, source, target, storePred)
	
	Abstract base class for source-target shortest path algorithms. 
	"""
	cdef Graph _G
	cdef vector[node] targets

	def __init__(self, *args, **namedargs):
		if type(self) == STSP:
			raise RuntimeError("Error, you may not use STSP directly, use a sub-class instead")

	def setSource(self, node newSource):
		"""
		setSource(newSource)

		Sets the source node.

		Parameters
		----------
		newSource : int
			The new source node.
		"""
		(<_STSP*>(self._this)).setSource(newSource)

	def setTarget(self, node newTarget):
		"""
		setTarget(newTarget)

		Sets the target node.

		Parameters
		----------
		newTarget : int
			The new target node.
		"""
		(<_STSP*>(self._this)).setTarget(newTarget)

	def setTargets(self, vector[node] newTargets):
		"""
		setTargets(newTargets)

		Sets multiple target nodes.

		Parameters
		----------
		newTargets : list(int)
			The new target nodes.
		"""
		self.targets = newTargets
		(<_STSP*>(self._this)).setTargets(newTargets.begin(), newTargets.end())

	def getPath(self):
		"""
		getPath()
		
		Returns a shortest path from the source node to the target node (without
		including them). Note: the shortest path can be constructed only if the
		algorithm is executed with @a storePred set to true.

		Returns
		-------
		list(int)
			A shortest path from the source node to the target node.
		"""
		return (<_STSP*>(self._this)).getPath()

	def getPredecessors(self):
		"""
		getPredecessors()

		Returns the predecessor nodes from the target node to the source node,
		Note: predecessors are stored only if the algorithm is executed with
		storePred set to true.

		Returns
		-------
		list(int)
			The list of predecessors from a target to a source.
		"""
		return (<_STSP*>(self._this)).getPredecessors()

	def getDistance(self):
		"""
		getDistance()

		Returns the distance from the source node to the target node

		Returns
		-------
		float
			The distance from source to the target node.
		"""
		return (<_STSP*>(self._this)).getDistance()

	def getDistances(self):
		"""
		getDistances()

		In case of multiple target nodes: returns the distance from the source node to the target
		nodes.

		Returns
		-------
		list(float)
			Distances from the source to the target nodes.
		"""
		tmap = (<_STSP*>(self._this)).getTargetIndexMap()
		distances = (<_STSP*>(self._this)).getDistances()
		return [distances[tmap[u]] for u in self.targets]

cdef extern from "<networkit/distance/SSSP.hpp>":

	cdef cppclass _SSSP "NetworKit::SSSP"(_Algorithm):
		_SSSP(_Graph G, node source, bool_t storePaths, bool_t storeNodesSortedByDistance, node target) except +
		vector[edgeweight] &getDistances() except +
		edgeweight distance(node t) except +
		vector[node] getPredecessors(node t) except +
		vector[node] getPath(node t, bool_t forward) except +
		set[vector[node]] getPaths(node t, bool_t forward) except +
		vector[node] getNodesSortedByDistance() except +
		double _numberOfPaths(node t) except +
		void setSource(node s) except +
		void setTarget(node t) except +

cdef class SSSP(Algorithm):
	""" 
	SSSP(G, source, storePaths, storeNodesSortedByDistance, target)

	Base class for single source shortest path algorithms. 
	"""
	cdef Graph _G

	def __init__(self, *args, **namedargs):
		if type(self) == SSSP:
			raise RuntimeError("Error, you may not use SSSP directly, use a sub-class instead")

	def getDistances(self, asarray=None):
		"""
		getDistances(asarray=None)

		Returns a list of weighted distances from the source node, i.e. the
 	 	length of the shortest path from the source node to any other node.

		Parameters
		----------
		asarray : optional
			Return the result as a numpy array. Default: None

 	 	Returns
 	 	-------
		list or np.ndarray
 	 		The weighted distances from the source node to any other node in the graph.
		"""
		return maybe_asarray_1d(&(<_SSSP*>(self._this)).getDistances(), asarray)

	def distance(self, t):
		"""
		distance(t)

		Returns the distance from the source node to t.

		Parameters
		----------
		t : int
			Target node.

		Returns
		-------
		float
			Distance from the source node to t.
		"""
		return (<_SSSP*>(self._this)).distance(t)

	def getPredecessors(self, t):
		"""
		getPredecessors(t)

		Returns the predecessor nodes of t on all shortest paths from source
		to t.

		Parameters
		----------
		t : int
			Target node.

		Returns
		-------
		list
			The predecessors of t on all shortest paths from source to t.
		"""
		return (<_SSSP*>(self._this)).getPredecessors(t)

	def getPath(self, t, forward=True):
		"""
		getPath(t, forward=True)

		Returns a shortest path from source to t and an empty path if source and t
		are not connected.

		Parameters
		----------
		t : int
			Target node.
		forward : bool, optional
			If True (default) the path is directed from source to t, otherwise the path
			is reversed.

		Returns
		-------
		list
			A shortest path from source to t or an empty path.
		"""
		return (<_SSSP*>(self._this)).getPath(t, forward)

	def getPaths(self, t, forward=True):
		"""
		getPaths(t, forward=True)

		Returns all shortest paths from source to t and an empty set if source
		and t are not connected.

		Parameters
		----------
		t : int
			Target node.
		forward : bool, optional
			If True (default) the path is directed from source to
			t, otherwise the path is reversed.

		Returns
		-------
			All shortest paths from source node to target node t.
		"""
		cdef set[vector[node]] paths = (<_SSSP*>(self._this)).getPaths(t, forward)
		result = []
		for elem in paths:
			result.append(list(elem))
		return result

	def getNodesSortedByDistance(self):
		""" 
		getNodesSortedByDistance()
		
		Returns a list of nodes ordered in increasing distance from the source.

		For this functionality to be available, storeNodesSortedByDistance has to be set to true in the constructor.
		There are no guarantees regarding the ordering of two nodes with the same distance to the source.

		Returns
		-------
		list
			Nodes ordered in increasing distance from the source.
		"""
		return (<_SSSP*>(self._this)).getNodesSortedByDistance()

	def numberOfPaths(self, t):
		"""
		numberOfPaths(t)

		Returns the number of paths from the source node to t.

		Parameters
		----------
		t : int
			Target node.

		Returns
		-------
		int
			The number of paths from the source node to t.
		"""
		return (<_SSSP*>(self._this))._numberOfPaths(t)

	def setSource(self, s not None):
		"""
		setSource(s)

		Sets a new source node.

		Parameters
		----------
		s : int
			New source node.
		"""
		(<_SSSP*>(self._this)).setSource(s)

	def setTarget(self, t not None):
		"""
		setTarget(t)

		Sets a new target node.

		Parameters
		----------
		t : int
			New target node.
		"""
		(<_SSSP*>(self._this)).setTarget(t)

cdef extern from "<networkit/distance/DynSSSP.hpp>":

	cdef cppclass _DynSSSP "NetworKit::DynSSSP"(_SSSP, _DynAlgorithm):
		_DynSSSP(_Graph G, node source, bool_t storePaths, bool_t storeStack, node target) except +
		bool_t modified() except +
		void setTargetNode(node t) except +

cdef class DynSSSP(SSSP, DynAlgorithm):
	""" 
	DynSSSP(G, source, storePredecessors, target)

	Base class for single source shortest path algorithms in dynamic graphs. 
	"""
	def __init__(self, *args, **namedargs):
		if type(self) == SSSP:
			raise RuntimeError("Error, you may not use DynSSSP directly, use a sub-class instead")

	def modified(self):
		""" 
		modified()

		Returns True or False depending on whether the node previoulsy specified with 
		setTargetNode(t) has been modified by the update or not.

		Returns
		-------
		bool
			Indicator for whether the target node was modified or not.
		"""
		return (<_DynSSSP*>(self._this)).modified()

	def setTargetNode(self, t):
		""" 
		setTargetNode(t)

		Set a target node to be observed during the update. If a node t is set as
		target before the update, the function modified() will return True or False
		depending on whether node t has been modified by the update.

		Parameters
		----------
		t : int
			Target node to be observed during update.
		"""
		(<_DynSSSP*>(self._this)).setTargetNode(t)

cdef extern from "<networkit/distance/AdamicAdarDistance.hpp>":

	cdef cppclass _AdamicAdarDistance "NetworKit::AdamicAdarDistance":
		_AdamicAdarDistance(const _Graph& G) except +
		void preprocess() except +
		double distance(node u, node v) except +
		vector[double] &getEdgeScores() except +

cdef class AdamicAdarDistance:
	"""
	AdamicAdarDistance(G)

	Calculate the adamic adar similarity.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""
	cdef _AdamicAdarDistance* _this
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _AdamicAdarDistance(G._this)

	def __dealloc__(self):
		del self._this

	def preprocess(self):
		self._this.preprocess()

	def getAttribute(self):
		"""
		getAttribute()

		Get the Adamic Adar similiraty score for every edge.

		Returns
		-------
		list(float)
			Adamic Adar similiraty score for every edge.

		"""
		#### TODO: convert distance to similarity!?! ####
		return self._this.getEdgeScores()

	def distance(self, node u, node v):
		"""
		distance(self, u, v)

		Calculate the distance from node u to node v.

		Parameters
		----------
		u : int
			Source node
		v : int
			Target node

		Returns
		-------
		float
			Distance from node u to node v.
		"""
		return self._this.distance(u, v)

cdef extern from "<networkit/distance/Diameter.hpp>" namespace "NetworKit":

	cdef enum _DiameterAlgo "NetworKit::DiameterAlgo":
		AUTOMATIC,
		EXACT,
		ESTIMATED_RANGE,
		ESTIMATED_SAMPLES,
		ESTIMATED_PEDANTIC,

class DiameterAlgo(object):
	AUTOMATIC = _DiameterAlgo.AUTOMATIC
	EXACT = _DiameterAlgo.EXACT
	ESTIMATED_RANGE = _DiameterAlgo.ESTIMATED_RANGE
	ESTIMATED_SAMPLES = _DiameterAlgo.ESTIMATED_SAMPLES
	ESTIMATED_PEDANTIC = _DiameterAlgo.ESTIMATED_PEDANTIC
	Automatic = AUTOMATIC # this + following added for backwards compatibility
	Exact = EXACT
	EstimatedRange = ESTIMATED_RANGE
	EstimatedSamples = ESTIMATED_SAMPLES
	EstimatedPedantic = ESTIMATED_PEDANTIC


cdef extern from "<networkit/distance/Diameter.hpp>" namespace "NetworKit::Diameter":

	cdef cppclass _Diameter "NetworKit::Diameter"(_Algorithm):
		_Diameter(_Graph G, _DiameterAlgo algo, double error, count nSamples) except +
		pair[count, count] getDiameter() except + nogil

cdef class Diameter(Algorithm):
	"""
	Diameter(G, algo = networkit.DiameterAlgo.AUTOMATIC, error = -1., nSamples = 0)

	Calculate the Diameter of the graph based different possible algorithms.

	Parameter :code:`algo` can be one of the following:

	- networkit.distance.DiameterAlgo.AUTOMATIC
	- networkit.distance.DiameterAlgo.EXACT
	- networkit.distance.DiameterAlgo.ESTIMATED_RANGE
	- networkit.distance.DiameterAlgo.ESTIMATED_SAMPLES
	- networkit.distance.DiameterAlgo.ESTIMATED_PEDANTIC

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	algo : networkit.distance.DiameterAlgo, optional
		Algorithm which should be used for diameter computation. Default: networkit.distance.DiameterAlgo.AUTOMATIC
	error : float, optional
		Possible error used for diameter algorithm EstimatedRange. Default: -1
	nSamples : int, optional
		Number of samples (influencing the quality of the output) used for diameter algorithm EstimatedSamples. Default: 0
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, algo = DiameterAlgo.AUTOMATIC, error = -1., nSamples = 0):
		self._G = G
		self._this = new _Diameter(G._this, algo, error, nSamples)

	def getDiameter(self):
		"""
		getDiameter()

		Returns
		-------
		float
			Diameter of the graph.
		"""
		return (<_Diameter*>(self._this)).getDiameter()

cdef extern from "<networkit/distance/Eccentricity.hpp>" namespace "NetworKit::Eccentricity":

	pair[node, count] getValue(_Graph G, node v) except +

cdef class Eccentricity:
	"""
	Eccentricity()
	
	The eccentricity of a node u is defined as the distance to the farthest node from node u. In other words, it is the longest shortest-path starting from node u.
	"""

	@staticmethod
	def getValue(Graph G, v):
		"""
		getValue(G, v)

		Get eccentricity value of node v from graph G.

		Parameters
		----------
		G : networkit.Graph
			The input graph.		

		Returns
		-------
		tuple(int, float)
			First index is the farthest node v from u, and the second index is the length of the shortest path from u to v.
		"""
		return getValue(G._this, v)

cdef extern from "<networkit/distance/EffectiveDiameterApproximation.hpp>" namespace "NetworKit::EffectiveDiameterApproximation":

	cdef cppclass _EffectiveDiameterApproximation "NetworKit::EffectiveDiameterApproximation"(_Algorithm):
		_EffectiveDiameterApproximation(_Graph& G, double ratio, count k, count r) except +
		double getEffectiveDiameter() except +

cdef class EffectiveDiameterApproximation(Algorithm):
	"""
	EffectiveDiameterApproximation(G, ratio=0.9, k=64, r=7)

	Calculates the effective diameter of a graph.
	The effective diameter is defined as the number of edges on average to reach a given ratio of all other nodes.

	Implementation after the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining in Massive Graphs"[1]

	[1] by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf

	Parameters
	----------
	G : networkit.Graph
		The graph.
	ratio : float, optional
		The percentage of nodes that shall be within stepwidth, default = 0.9
	k : int, optional
		Number of parallel approximations, bigger k -> longer runtime, more precise result; default = 64
	r : int, optional
		Number of additional bits, important in tiny graphs; default = 7
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, double ratio=0.9, count k=64, count r=7):
		self._G = G
		self._this = new _EffectiveDiameterApproximation(G._this, ratio, k, r)

	def getEffectiveDiameter(self):
		"""
		getEffectiveDiameter()

		Returns
		-------
		float
			The approximated effective diameter
		"""
		return (<_EffectiveDiameterApproximation*>(self._this)).getEffectiveDiameter()

cdef extern from "<networkit/distance/EffectiveDiameter.hpp>" namespace "NetworKit::EffectiveDiameter":

	cdef cppclass _EffectiveDiameter "NetworKit::EffectiveDiameter"(_Algorithm):
		_EffectiveDiameter(_Graph& G, double ratio) except +
		double getEffectiveDiameter() except +

cdef class EffectiveDiameter(Algorithm):
	"""
	EffectiveDiameter(G, ratio=0.9)

	Calculates the effective diameter of a graph.
	The effective diameter is defined as the number of edges on average to reach a given ratio of all other nodes.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	ratio : float, optional
		The percentage of nodes that shall be within stepwidth; default = 0.9
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, double ratio=0.9):
		self._G = G
		self._this = new _EffectiveDiameter(G._this, ratio)

	def getEffectiveDiameter(self):
		"""
		getEffectiveDiameter()

		Returns
		-------
		float
			The effective diameter
		"""
		return (<_EffectiveDiameter*>(self._this)).getEffectiveDiameter()

cdef extern from "<networkit/distance/HopPlotApproximation.hpp>" namespace "NetworKit::HopPlotApproximation":

	cdef cppclass _HopPlotApproximation "NetworKit::HopPlotApproximation"(_Algorithm):
		_HopPlotApproximation(_Graph& G, count maxDistance, count k, count r) except +
		map[count, double] &getHopPlot() except +

cdef class HopPlotApproximation(Algorithm):
	"""
	HopPlotApproximation(G, maxDistance=0, k=64, r=7)

	Computes an approxmation of the hop-plot of a given graph.
	The hop-plot is the set of pairs (d, g(g)) for each natural number d
	and where g(d) is the fraction of connected node pairs whose shortest connecting path has length at most d.

	Implementation after the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining in Massive Graphs"[1]

	[1] by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf

	Parameters
	----------
	G : networkit.Graph
		The graph.
	maxDistance : float, optional
		Maximum distance between considered nodes set to 0 or negative to get the hop-plot 
		for the entire graph so that each node can reach each other node.
	k : int, optional
		Number of parallel approximations, bigger k -> longer runtime, more precise result; default = 64
	r : int, optional
		Number of additional bits, important in tiny graphs; default = 7
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, count maxDistance=0, count k=64, count r=7):
		self._G = G
		self._this = new _HopPlotApproximation(G._this, maxDistance, k, r)

	def getHopPlot(self):
		"""
		getHopPlot()

		Returns the approximated hop-plot of the graph.

		Returns
		-------
		dict(int ``:`` float)
			Number of connected nodes for each distance
		"""
		cdef map[count, double] hp = (<_HopPlotApproximation*>(self._this)).getHopPlot()
		result = dict()
		for elem in hp:
			result[elem.first] = elem.second
		return result

cdef extern from "<networkit/distance/NeighborhoodFunction.hpp>" namespace "NetworKit::NeighborhoodFunction":

	cdef cppclass _NeighborhoodFunction "NetworKit::NeighborhoodFunction"(_Algorithm):
		_NeighborhoodFunction(_Graph& G) except +
		vector[count] &getNeighborhoodFunction() except +

cdef class NeighborhoodFunction(Algorithm):
	"""
	NeighborhoodFunction(G)

	Computes the neighborhood function exactly.
	The neighborhood function N of a graph G for a given distance t is defined
	as the number of node pairs (u,v) that can be reached within distance t.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None):
		self._G = G
		self._this = new _NeighborhoodFunction(G._this)

	def getNeighborhoodFunction(self):
		"""
		getNeighborhoodFunction()

		Returns the neighborhood function of the graph.

		Returns
		-------
		list(int)
			The i-th element denotes the number of node pairs that have a distance at most (i+1).
		"""
		return (<_NeighborhoodFunction*>(self._this)).getNeighborhoodFunction()

cdef extern from "<networkit/distance/NeighborhoodFunctionApproximation.hpp>" namespace "NetworKit::NeighborhoodFunctionApproximation":

	cdef cppclass _NeighborhoodFunctionApproximation "NetworKit::NeighborhoodFunctionApproximation"(_Algorithm):
		_NeighborhoodFunctionApproximation(_Graph& G, count k, count r) except +
		vector[count] &getNeighborhoodFunction() except +

cdef class NeighborhoodFunctionApproximation(Algorithm):
	"""
	NeighborhoodFunctionApproximation(G, k=64, r=7)

	Computes an approximation of the neighborhood function.
	The neighborhood function N of a graph G for a given distance t is defined
	as the number of node pairs (u,v) that can be reached within distance t.

	Implementation after the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining in Massive Graphs"[1]

	[1] by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf

	Parameters
	----------
	G : networkit.Graph
		The graph.
	k : int, optional
		Number of approximations, bigger k -> longer runtime, more precise result; default = 64
	r : int, optional
		Number of additional bits, important in tiny graphs; default = 7
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, count k=64, count r=7):
		self._G = G
		self._this = new _NeighborhoodFunctionApproximation(G._this, k, r)

	def getNeighborhoodFunction(self):
		"""
		getNeighborhoodFunction()

		Returns the neighborhood function of the graph.

		Returns
		-------
		list(int)
			The i-th element denotes the number of node pairs that have a distance at most (i+1).
		"""
		return (<_NeighborhoodFunctionApproximation*>(self._this)).getNeighborhoodFunction()

cdef extern from "<networkit/distance/Volume.hpp>" namespace "NetworKit::Volume":

	double volume(const _Graph G, const double r, const count samples) except + nogil
	vector[double] volume(const _Graph G, const vector[double] r, const count samples) except + nogil

cdef class Volume:
	"""
	Volume()
	
	The volume of a graph and its meaning is explained in the following publication:

	Franz-Benjamin Mocnik: "The Polynomial Volume Law of Complex Networks in
	the Context of Local and Global Optimization", Scientific Reports 8(11274)
	2018. doi: 10.1038/s41598-018-29131-0
	"""

	@staticmethod
	def volume(Graph G, r, count samples=500):
		"""
		volume(G, r, samples=500)

		Number of nodes within a given radius. If the radius is a list containing many radii,
		a list containing the number for every radius is returned.

		Parameters
		----------
		G : networkit.Graph
			the graph
		r : float
			the radius (or radii)
		samples : int, optional
			The number of samples. Default: 500

		Returns
		-------
		float
			Number of nodes within a given radius.
		list(float)
			Number of nodes within every given radii.
		"""
		cdef double _r
		cdef vector[double] _rs
		cdef double _v
		cdef vector[double] _vs
		def is_number(s):
			try:
				float(s)
				return True
			except ValueError:
				return False
		if type(r) is float or type(r) is int:
			_r = r
			with nogil:
				_v = volume(<_Graph> G._this, <double> _r, <count> samples)
			return _v
		elif type(r) is list and all(is_number(item) for item in r):
			_rs = r
			with nogil:
				_vs = volume(<_Graph> G._this, <vector[double]> _rs, <count> samples)
			return _vs
		else:
			pass

cdef extern from "<networkit/distance/JaccardDistance.hpp>":

	cdef cppclass _JaccardDistance "NetworKit::JaccardDistance":
		_JaccardDistance(const _Graph& G, const vector[count]& _triangles) except +
		void preprocess() except +
		vector[double] &getEdgeScores() except +

cdef class JaccardDistance:
	"""
	JaccardDistance(G, triangles)

	The Jaccard distance measure assigns to each edge the jaccard coefficient
	of the neighborhoods of the two adjacent nodes.

	Parameters:
	-----------
	G : networkit.Graph
		The graph to calculate Jaccard distances for.
	triangles : list(int)
		Previously calculated edge triangle counts.
	"""

	cdef _JaccardDistance* _this
	cdef Graph _G
	cdef vector[count] _triangles

	def __cinit__(self, Graph G, vector[count] triangles):
		self._G = G
		self._triangles = triangles
		self._this = new _JaccardDistance(G._this, self._triangles)

	def __dealloc__(self):
		del self._this

	def getAttribute(self):
		"""
		getAttribute()

		Get the Jaccard distance for every edge.

		Returns
		-------
		list(float)
			Jaccard distance for every edge.

		"""
		return self._this.getEdgeScores()

cdef class JaccardSimilarityAttributizer:
	"""
	JaccardSimilarityAtrributizer(G, triangles)

	The Jaccard similarity measure assigns to each edge (1 - the jaccard coefficient
	of the neighborhoods of the two adjacent nodes).

	Parameters
	----------
	G : networkit.Graph
		The graph to calculate Jaccard similarities for.
	triangles : list(int)
		Previously calculated edge triangle counts.
	"""

	cdef _JaccardDistance* _this
	cdef Graph _G
	cdef vector[count] _triangles

	def __cinit__(self, Graph G, vector[count] triangles):
		self._G = G
		self._triangles = triangles
		self._this = new _JaccardDistance(G._this, self._triangles)

	def __dealloc__(self):
		del self._this

	def getAttribute(self):
		"""
		getAttribute()

		Get the Jaccard similiraty score for every edge.

		Returns
		-------
		list(float)
			Jaccard similiraty score for every edge.

		"""
		#convert distance to similarity
		self._this.preprocess()
		return [1 - x for x in self._this.getEdgeScores()]

cdef extern from "<networkit/distance/AlgebraicDistance.hpp>":

	cdef cppclass _AlgebraicDistance "NetworKit::AlgebraicDistance":
		_AlgebraicDistance(_Graph G, count numberSystems, count numberIterations, double omega, index norm, bool_t withEdgeScores) except +
		void preprocess() except +
		double distance(node, node) except +
		vector[double] &getEdgeScores() except +


cdef class AlgebraicDistance:
	"""
	AlgebraicDistance(G, numberSystems=10, numberIterations=30, omega=0.5, norm=0, withEdgeScores=False)

	Algebraic distance assigns a distance value to pairs of nodes
	according to their structural closeness in the graph.
	Algebraic distances will become small within dense subgraphs.

	Parameters
	----------
	G : networkit.Graph
		The graph to calculate Jaccard distances for.
	numberSystems : int, optional
	 	Number of vectors/systems used for algebraic iteration. Default: 10
	numberIterations : int, optional
	 	Number of iterations in each system. Default: 30
	omega : float, optional
	 	Attenuation factor in [0,1] influencing convergence speed. Default: 0.5
	norm : int, optional
		The norm factor of the extended algebraic distance. Default: 0
	withEdgeScores : bool, optional
		Calculate array of scores for edges {u,v} that equal ad(u,v). Default: False
	"""

	cdef _AlgebraicDistance* _this
	cdef Graph _G

	def __cinit__(self, Graph G, count numberSystems=10, count numberIterations=30, double omega=0.5, index norm=0, bool_t withEdgeScores=False):
		self._G = G
		self._this = new _AlgebraicDistance(G._this, numberSystems, numberIterations, omega, norm, withEdgeScores)

	def __dealloc__(self):
		del self._this

	def preprocess(self):
		self._this.preprocess()
		return self

	def distance(self, node u, node v):
		return self._this.distance(u, v)

	def getEdgeScores(self):
		return self._this.getEdgeScores()

cdef extern from "<networkit/distance/CommuteTimeDistance.hpp>":

	cdef cppclass _CommuteTimeDistance "NetworKit::CommuteTimeDistance"(_Algorithm):
		_CommuteTimeDistance(_Graph G, double tol) except +
		void runApproximation() except +
		void runParallelApproximation() except +
		double distance(node, node) except +
		double runSinglePair(node, node) except +
		double runSingleSource(node) except +


cdef class CommuteTimeDistance(Algorithm):
	""" 
	CommuteTimeDistance(G, tol=0.1)
	
	Computes the Euclidean Commute Time Distance (ECTD) between each pair of nodes for an undirected unweighted graph.

	CommuteTimeDistance(G)

	Parameters
	----------
	G : networkit.Graph
		The graph.
	tol: float, optional
		Tolerance for computation (higher tolerance leads to faster running times). Default: 0.1
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G, double tol = 0.1):
		self._G = G
		self._this = new _CommuteTimeDistance(G._this, tol)

	def runApproximation(self):
		""" 
		runApproximation()
		
		Computes approximation of the ECTD. 
		"""
		return (<_CommuteTimeDistance*>(self._this)).runApproximation()

	def runParallelApproximation(self):
		""" 
		runParallelApproximation()
		
		Computes approximation (in parallel) of the ECTD. 
		"""
		return (<_CommuteTimeDistance*>(self._this)).runParallelApproximation()

	def distance(self, u, v):
		"""
		distance(u, v)

		Returns the ECTD between node u and node v.

		Parameters
		----------
		u : int
			Index of node u.
		v : int
			Index of node v.

		Returns
		-------
		float
			ECTD between node u and v.
		"""
		return (<_CommuteTimeDistance*>(self._this)).distance(u, v)

	def runSinglePair(self, u, v):
		"""
		runSinglePair(u, v)
		
		Returns the ECTD between node u and node v, without preprocessing.

		Parameters
		----------
		u : int
			Index of node u.
		v : int
			Index of node v.

		Returns
		-------
		float
			ECTD (without preprocessing) between node u and v.
		"""
		return (<_CommuteTimeDistance*>(self._this)).runSinglePair(u, v)

	def runSingleSource(self, u):
		"""
		runSingleSource(u)
		
		Returns the sum of the ECTDs from u, without preprocessing.

		Parameters
		----------
		u : int
			Index of node u.

		Returns
		-------
		float
			Sum of the ECTDs from u, without preprocessing.
		"""
		return (<_CommuteTimeDistance*>(self._this)).runSingleSource(u)

cdef extern from "<networkit/distance/NeighborhoodFunctionHeuristic.hpp>" namespace "NetworKit::NeighborhoodFunctionHeuristic::SelectionStrategy":

	enum _SelectionStrategy "NetworKit::NeighborhoodFunctionHeuristic::SelectionStrategy":
		RANDOM = 0
		SPLIT = 1

class SelectionStrategy(object):
	RANDOM = _SelectionStrategy.RANDOM
	SPLIT = _SelectionStrategy.SPLIT

cdef extern from "<networkit/distance/NeighborhoodFunctionHeuristic.hpp>" namespace "NetworKit::NeighborhoodFunctionHeuristic":

	cdef cppclass _NeighborhoodFunctionHeuristic "NetworKit::NeighborhoodFunctionHeuristic"(_Algorithm):
		_NeighborhoodFunctionHeuristic(_Graph& G, const count nSamples, const _SelectionStrategy strategy) except +
		vector[count] &getNeighborhoodFunction() except +

cdef class NeighborhoodFunctionHeuristic(Algorithm):
	"""
	NeighborhoodFunctionHeuristic(G, nSamples=0, strategy=SelectionStrategy.SPLIT)

	Computes a heuristic of the neighborhood function.
	The algorithm runs nSamples breadth-first searches and scales the results up to the actual amount of nodes.

	Parameter :code:`strategy` can be one of the following:

	- networkit.distance.SelectionStrategy.RANDOM
	- networkit.distance.SelectionStrategy.SPLIT

	Parameters
	----------
	G : networkit.Graph
		The graph.
	nSamples : int, optional
		The amount of samples, set to zero for heuristic of max(sqrt(m), 0.15*n). Default: 0
	strategy : networkit.distance.SelectionStrategy, optional
		The strategy to select the samples, accepts RANDOM (0) or SPLIT (1). Default: networkit.distance.SelectionStrategy.SPLIT
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, count nSamples=0, strategy=SelectionStrategy.SPLIT):
		self._G = G
		self._this = new _NeighborhoodFunctionHeuristic(G._this, nSamples, strategy)

	def getNeighborhoodFunction(self):
		"""
		getNeighborhoodFunction()

		Returns the neighborhood function of the graph.

		Returns
		-------
		list(int)
			The i-th element denotes the number of node pairs that have a distance at most (i+1).
		"""
		return (<_NeighborhoodFunctionHeuristic*>(self._this)).getNeighborhoodFunction()

cdef extern from "<networkit/distance/APSP.hpp>":

	cdef cppclass _APSP "NetworKit::APSP"(_Algorithm):
		_APSP(_Graph G) except +
		vector[vector[edgeweight]] &getDistances() except +
		edgeweight getDistance(node u, node v) except +

cdef class APSP(Algorithm):
	""" 
	APSP(G)

	All-Pairs Shortest-Paths algorithm (implemented running Dijkstra's algorithm from each node, or BFS if G is unweighted).
	Computes all pairwise shortest-path distances in G.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _APSP(G._this)

	def __dealloc__(self):
		self._G = None

	def getDistances(self, asarray=None):
		"""
		getDistances(asarray=None)

		Returns a vector of vectors of distances between each node pair.

		Parameters
		----------
		asarray : optional
			Return the result as a numpy array. Default: Falsy.

		Returns
		-------
		list(list(float)) or np.ndarray
			The shortest-path distances from each node to any other node in the graph.
		"""
		return maybe_asarray_2d(&(<_APSP*>(self._this)).getDistances(), asarray)

	def getDistance(self, node u, node v):
		"""
		getDistance(u, v)
		
		Returns the length of the shortest path from source u to target v.

		Parameters
		----------
		u : node
			Index of source node u.
		v : node
			Index of target node v.

		Returns
		-------
		int or float
			The distance from u to v. Returned value is of type int, if the graph is unweighted - otherwise the return
			type is float.
		"""
		return (<_APSP*>(self._this)).getDistance(u, v)

cdef extern from "<networkit/distance/SPSP.hpp>":

	cdef cppclass _SPSP "NetworKit::SPSP"(_Algorithm):
		_SPSP(_Graph G, vector[node].iterator sourcesFirst, vector[node].iterator sourcesLast) except +
		_SPSP(_Graph G, vector[node].iterator sourcesFirst, vector[node].iterator sourcesLast,
			vector[node].iterator targetsFirst, vector[node].iterator targetsLast) except +
		vector[vector[edgeweight]] &getDistances() except +
		edgeweight getDistance(node u, node v) except +
		void setSources(vector[node].iterator sourcesFirst, vector[node].iterator sourcesLast)
		void setTargets(vector[node].iterator targetsFirst, vector[node].iterator targetsLast)

cdef class SPSP(Algorithm):
	"""
	SPSP(G, sources)
	
	Some-Pairs Shortest-Paths algorithm (implemented running Dijkstra's algorithm from each source
	node, or BFS if G is unweighted).
	Computes pairwise shortest-path distances from the source nodes to all the nodes in G.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	sources : list(int)
		Set of source nodes.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, vector[node] sources, vector[node] targets = []):
		self._G = G
		if not targets.empty():
			self._this = new _SPSP(G._this, sources.begin(), sources.end(), targets.begin(), targets.end())
		else:
			self._this = new _SPSP(G._this, sources.begin(), sources.end())

	def __dealloc__(self):
		self._G = None

	def getDistances(self, asarray=None):
		"""
		getDistances(asarray=None)

		Returns a list of lists of distances between each source node
		and either all the other nodes (if no targets have been specified) or all target nodes.

		Parameters
		----------
		asarray : optional
			Return the result as a numpy array. Default: Falsy.

		Returns
		-------
		list(list(float)) or np.ndarray
			The shortest-path distances from each source node to the target nodes (if target nodes
			have been specified), or any other node in the graph.

		Returns
		-------
		list(list(float)) or np.ndarray
			The shortest-path distances from each source node to any other node
			in the graph.
		"""
		return maybe_asarray_2d(&(<_SPSP*>self._this).getDistances(), asarray)

	def getDistance(self, node u, node v):
		"""
		getDistance(u, v)

		Returns the length of the shortest path from source u to target v.

		Parameters
		----------
		u : node
			Index of source node.
		v : node
			Index of target node.

		Returns
		-------
		float
			The distance from u to v.
		"""
		return (<_SPSP*>self._this).getDistance(u, v)

	def setSources(self, vector[node] sources):
		""" 
		setSources(sources)
		
		Sets the source nodes.

		Parameters
		----------
		sources : list(int)
			List of the new source nodes.
		"""
		(<_SPSP*>self._this).setSources(sources.begin(), sources.end())

	def setTargets(self, vector[node] targets):
		"""
		setTargets(targets)
		
		Sets the target nodes.

		Parameters
		----------
		targets : list(int)
			List of the new target nodes.
		"""
		(<_SPSP*>self._this).setTargets(targets.begin(), targets.end())

cdef extern from "<networkit/distance/DynAPSP.hpp>":

	cdef cppclass _DynAPSP "NetworKit::DynAPSP"(_APSP, _DynAlgorithm):
		_DynAPSP(_Graph G) except +

cdef class DynAPSP(APSP, DynAlgorithm):
	""" 
	DynAPSP(G)
	
	All-Pairs Shortest-Paths algorithm for dynamic graphs.
	Computes all pairwise shortest-path distances in G.

	Parameters
	----------
	G : networkit.Graph
		The graph.
		"""
	def __init__(self, Graph G):
		self._G = G
		self._this = new _DynAPSP(G._this)



cdef extern from "<networkit/distance/BFS.hpp>":

	cdef cppclass _BFS "NetworKit::BFS"(_SSSP):
		_BFS(_Graph G, node source, bool_t storePaths, bool_t storeNodesSortedByDistance, node target) except +

cdef class BFS(SSSP):
	""" 
	BFS(G, source, storePaths=True, storeNodesSortedByDistance=False, target=None)
	
	Simple breadth-first search on a Graph from a given source.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	source : int
		The source node of the breadth-first search.
	storePaths : bool, optional
		Controls whether to store paths and number of paths. Default: True
	storeNodesSortedByDistance : bool, optional
		Controls whether to store nodes sorted by distance. Default: False
	target: int or None, optional
		Terminate search when the target has been reached. In default-mode, this target is set to None.
	"""

	def __cinit__(self, Graph G, source, storePaths=True, storeNodesSortedByDistance=False, target=none):
		self._G = G
		self._this = new _BFS(G._this, source, storePaths, storeNodesSortedByDistance, target)

cdef extern from "<networkit/distance/Dijkstra.hpp>":

	cdef cppclass _Dijkstra "NetworKit::Dijkstra"(_SSSP):
		_Dijkstra(_Graph G, node source, bool_t storePaths, bool_t storeNodesSortedByDistance, node target) except +

cdef class Dijkstra(SSSP):
	""" 
	Dijkstra(G, source, storePaths=True, storeNodesSortedByDistance=False, target=None)
	
	Dijkstra's SSSP algorithm. Returns list of weighted distances from node source, i.e. the length of the shortest path from source to
	any other node.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	source : int
		The source node of the Dijkstra search.
	storePaths : bool, optional
		Controls whether to store paths and number of paths. Default: True
	storeNodesSortedByDistance : bool, optional
		Controls whether to store nodes sorted by distance. Default: False
	target: int or None, optional
		Terminate search when the target has been reached. In default-mode, this target is set to None.
	"""
	def __cinit__(self, Graph G, source, storePaths=True, storeNodesSortedByDistance=False, node target=none):
		self._G = G
		self._this = new _Dijkstra(G._this, source, storePaths, storeNodesSortedByDistance, target)

cdef extern from "<networkit/distance/MultiTargetBFS.hpp>":
	cdef cppclass _MultiTargetBFS "NetworKit::MultiTargetBFS"(_STSP):
		_MultiTargetBFS(_Graph G, node source, vector[node].iterator targetsFirst, vector[node].iterator targetsLast) except +

cdef class MultiTargetBFS(STSP):
	"""
	MultiTargetBFS(G, source, targets)
	
	Simple breadth-first search on a Graph from a given source to multiple targets.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	source : int
		The source node of the breadth-first search.
	targets : list(int)
		List of target nodes.
	"""

	def __cinit__(self, Graph G, node source, vector[node] targets):
		self._G = G
		self._this = new _MultiTargetBFS(G._this, source, targets.begin(), targets.end())
		self.targets = targets

cdef extern from "<networkit/distance/MultiTargetDijkstra.hpp>":
	cdef cppclass _MultiTargetDijkstra "NetworKit::MultiTargetDijkstra"(_STSP):
		_MultiTargetDijkstra(_Graph G, node source, vector[node].iterator targetsFirst, vector[node].iterator targetsLast) except +

cdef class MultiTargetDijkstra(STSP):
	"""
	MultiTargetDijkstra(G, source, targets)
	
	Dijkstra's SSSP algorithm from a single source node to multiple target nodes.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	source : int
		The source node of the Dijkstra search.
	targets : list(int)
		List of target nodes.
	"""

	def __cinit__(self, Graph G, node source, vector[node] targets):
		self._G = G
		self._this = new _MultiTargetDijkstra(G._this, source, targets.begin(), targets.end())
		self.targets = targets

cdef extern from "<networkit/distance/DynBFS.hpp>":

	cdef cppclass _DynBFS "NetworKit::DynBFS"(_DynSSSP):
		_DynBFS(_Graph G, node source) except +

cdef class DynBFS(DynSSSP):
	""" 
	DynBFS(G, source)

	Dynamic version of BFS.	

	Parameters
	----------
	G : networkit.Graph
		The graph.
	source : int
		The source node of the breadth-first search.
	"""
	def __cinit__(self, Graph G, source):
		self._G = G
		self._this = new _DynBFS(G._this, source)

cdef extern from "<networkit/distance/DynDijkstra.hpp>":

	cdef cppclass _DynDijkstra "NetworKit::DynDijkstra"(_DynSSSP):
		_DynDijkstra(_Graph G, node source) except +

cdef class DynDijkstra(DynSSSP):
	""" 
	DynDijkstra(G, source)	
	
	Dynamic version of Dijkstra. Create DynDijkstra for G and a source node.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	source : int
		The source node of the Dijkstra search.

	"""
	def __cinit__(self, Graph G, source):
		self._G = G
		self._this = new _DynDijkstra(G._this, source)

cdef cppclass PathCallbackWrapper:
	void* callback
	__init__(object callback):
		this.callback = <void*>callback
	void cython_call_operator(vector[node] path):
		cdef bool_t error = False
		cdef string message
		try:
			(<object>callback)(path)
		except Exception as e:
			error = True
			message = stdstring("An Exception occurred, aborting execution of iterator: {0}".format(e))
		if (error):
			throw_runtime_error(message)

cdef extern from "<networkit/distance/BidirectionalBFS.hpp>":
	cdef cppclass _BidirectionalBFS "NetworKit::BidirectionalBFS"(_STSP):
		_BidirectionalBFS(_Graph G, node source, node target, bool_t storePred) except +

cdef class BidirectionalBFS(STSP):
	"""
	BidirectionalBFS(G, source, target, storePre=True)

	Implements a bidirectional breadth-first search on a graph from two given source and target nodes.
	Explores the graph from both the source and target nodes until the two explorations meet.

	Parameters
	----------

	G : networkit.Graph
		The input graph.
	source : int
		The source node.
	target : int
		The target node.
	storePred : bool, optional
		If True, the algorithm will also store the predecessors
		and reconstruct a shortest path from source and target. Default: True
	"""

	def __cinit__(self, Graph G, node source, node target, bool_t storePred=True):
		self._this = new _BidirectionalBFS(G._this, source, target, storePred)

cdef extern from "<networkit/distance/BidirectionalDijkstra.hpp>":
	cdef cppclass _BidirectionalDijkstra "NetworKit::BidirectionalDijkstra"(_STSP):
		_BidirectionalDijkstra(_Graph G, node source, node target, bool_t storePred) except +

cdef class BidirectionalDijkstra(STSP):
	"""
	BidirectionalDijkstra(G, source, target, storePred=True)

	Bidirectional implementation of the Dijkstra algorithm from
	two given source and target nodes.
	Explores the graph from both the source and target nodes until
	the two explorations meet.

	Parameters:
	-----------

	G : networkit.Graph
		The input graph.
	source : int
		The source node.
	target : int
		The target node.
	storePred : bool, optional
		If True, the algorithm will also store the predecessors
		and reconstruct a shortest path from source and target.
	"""

	def __cinit__(self, Graph G, node source, node target, bool_t storePred=True):
		self._this = new _BidirectionalDijkstra(G._this, source, target, storePred)

cdef extern from "<networkit/distance/AStar.hpp>":
	cdef cppclass _AStar "NetworKit::AStar"(_STSP):
		_AStar(_Graph G, vector[double] &heu, node source, node target, bool_t storePred) except +

cdef class AStar(STSP):
	"""
	AStar(G, heu, source, target, storePred=True)

	A* path-finding algorithm.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	heu : list(float)
		List of lower bounds of the distance of each node to the target.
	source : int
		The source node.
	target : int
		The target node.
	storePred : bool, optional
		If True, the algorithm will also store the predecessors
		and reconstruct a shortest path from source and target. Default: True
	"""

	cdef vector[double] heu
	def __cinit__(self, Graph G, vector[double] &heu, node source, node target, bool_t storePred=True):
		self.heu = heu
		self._this = new _AStar(G._this, self.heu, source, target, storePred)

cdef extern from "<networkit/reachability/AllSimplePaths.hpp>":

	cdef cppclass _AllSimplePaths "NetworKit::AllSimplePaths":
		_AllSimplePaths(_Graph G, node source, node target, count cutoff) except +
		void run() except + nogil
		count numberOfSimplePaths() except +
		vector[vector[node]] getAllSimplePaths() except +
		void forAllSimplePaths[Callback](Callback c) except +

cdef class AllSimplePaths:
	""" 
	AllSimplePaths(G, source, target, cutoff=None)

	Algorithm to compute all existing simple paths from a source node to a target node. The maximum length of the paths can be fixed through 'cutoff'.
	CAUTION: This algorithm could take a lot of time on large networks (many edges), especially if the cutoff value is high or not specified.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	source : int
		The source node.
	target : int
		The target node.
	cutoff : int or None, optional
		The maximum length of the simple paths. In default mode, there is no cutoff set.
	"""

	cdef _AllSimplePaths* _this
	cdef Graph _G

	def __cinit__(self,  Graph G, source, target, cutoff=none):
		self._G = G
		self._this = new _AllSimplePaths(G._this, source, target, cutoff)
		from warnings import warn
		warn("networkit.distance.AllSimplePaths is deprecated, use networkit.reachability.AllSimplePaths")

	def __dealloc__(self):
		del self._this

	def run(self):
		self._this.run()
		return self

	def numberOfSimplePaths(self):
		"""
		numberOfSimplePaths()

		Returns the number of simple paths.

		Returns
		-------
		int
			The number of simple paths.
		"""
		return self._this.numberOfSimplePaths()

	def getAllSimplePaths(self):
		"""
		getAllSimplePaths()

		Returns all the simple paths from source to target.

		Returns
		-------
		list(list(int))
			A list containing list of node indexes which represent all simple paths.
		"""
		return self._this.getAllSimplePaths()

	def forAllSimplePaths(self, object callback):
		""" 
		forAllSimplePaths(callback)
		
		More efficient path iterator. Iterates over all the simple paths.

		Parameters
		----------
		callback : object
			Any callable object that takes the parameter path
		"""
		cdef PathCallbackWrapper* wrapper
		try:
			wrapper = new PathCallbackWrapper(callback)
			self._this.forAllSimplePaths[PathCallbackWrapper](dereference(wrapper))
		finally:
			del wrapper

cdef extern from "<networkit/distance/ReverseBFS.hpp>":

	cdef cppclass _ReverseBFS "NetworKit::ReverseBFS"(_SSSP):
		_ReverseBFS(_Graph G, node source, bool_t storePaths, bool_t storeNodesSortedByDistance, node target) except +

cdef class ReverseBFS(SSSP):
	""" 
	ReverseBFS(G, source, storePaths=True, storeNodesSortedByDistance=False, target=None)

	Simple reverse breadth-first search on a Graph from a given source.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	source : int
		The source node of the breadth-first search.
	storePaths : bool, optional
		Controls whether to store paths and number of paths. Default: True
	storeNodesSortedByDistance : bool, optional
		Controls whether to store nodes sorted by distance. Default: False
	target: int or None, optional
		Terminate search when the target has been reached. In default-mode, this target is set to None.
	"""

	def __cinit__(self, Graph G, source, storePaths=True, storeNodesSortedByDistance=False, target=none):
		self._G = G
		self._this = new _ReverseBFS(G._this, source, storePaths, storeNodesSortedByDistance, target)

cdef extern from "<networkit/distance/PrunedLandmarkLabeling.hpp>":

	cdef cppclass _PrunedLandmarkLabeling "NetworKit::PrunedLandmarkLabeling"(_Algorithm):
		_PrunedLandmarkLabeling(_Graph G) except +
		count query(node u, node v) except +

cdef class PrunedLandmarkLabeling(Algorithm):
	"""
	PrunedLandmarkLabeling(G)

	Pruned Landmark Labeling algorithm based on the paper "Fast exact shortest-path distance
	queries on large networks by pruned landmark labeling" from Akiba et al., ACM SIGMOD 2013.
	The algorithm computes distance labels by performing pruned breadth-first searches from each
	vertex. Labels are used to quickly retrieve shortest-path distances between node pairs.
	Note: this algorithm only works for unweighted graphs.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _PrunedLandmarkLabeling(G._this)

	def __dealloc__(self):
		self._G = None

	def query(self, node u, node v):
		"""
		query(u, v)

		Returns the shortest-path distance between the two nodes.

		Parameters
		----------
		u : node
			Source node.
		v : node
			Target node.

		Returns
		-------
		int
			The shortest-path distances from the source node to the target node.
		"""
		return (<_PrunedLandmarkLabeling*>(self._this)).query(u, v)


cdef extern from "<networkit/distance/DynPrunedLandmarkLabeling.hpp>":

	cdef cppclass _DynPrunedLandmarkLabeling "NetworKit::DynPrunedLandmarkLabeling"(_Algorithm, _DynAlgorithm):
		_DynPrunedLandmarkLabeling(_Graph G) except +
		count query(node u, node v) except +

cdef class DynPrunedLandmarkLabeling(Algorithm, DynAlgorithm):
	"""
	DynPrunedLandmarkLabeling(G)

	Dynamic Pruned Landmark Labeling algorithm based on the paper "Fully
	Dynamic 2-Hop Cover Labeling " from D'Angelo et al., ACM JEA 2019. The
	algorithm computes distance labels by performing pruned breadth-first
	searches from each vertex. Distance labels can be updated efficiently
	after edge insertions.
	Note: this algorithm only works for unweighted graphs and only supports
	edge insertions.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _DynPrunedLandmarkLabeling(G._this)

	def __dealloc__(self):
		self._G = None

	def query(self, node u, node v):
		"""
		query(u, v)

		Returns the shortest-path distance between the two nodes.

		Parameters
		----------
		u : node
			Source node.
		v : node
			Target node.

		Returns
		-------
		int
			The shortest-path distances from the source node to the target node.
		"""
		return (<_DynPrunedLandmarkLabeling*>(self._this)).query(u, v)
