# distutils: language=c++

from cython.operator import dereference, preincrement

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.map cimport map
from libcpp cimport bool as bool_t
from libcpp.string cimport string
from libcpp.set cimport set

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node
ctypedef double edgeweight

from .base cimport _Algorithm, Algorithm
from .dynamics cimport _GraphEvent
from .graph cimport _Graph, Graph
from .helpers import stdstring

cdef extern from "<networkit/Globals.hpp>" namespace "NetworKit":

	index _none "NetworKit::none"

none = _none

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message)

cdef extern from "<networkit/distance/STSP.hpp>":
	cdef cppclass _STSP "NetworKit::STSP"(_Algorithm):
		_STSP(_Graph G, node source, node target, bool_t storePred) except +
		vector[node] getPath() except +
		vector[node] getPredecessors() except +
		edgeweight getDistance() except +

cdef class STSP(Algorithm):
	""" Abstract base class for source-target shortest path algorithms. """
	cdef Graph _G

	def __init__(self, *args, **namedargs):
		"""
		Creates the STSP class for a graph G, a source node, and a target node.

		Parameters:
		-----------
		G : networkit.Graph
			The graph.
		source : node
			The source node.
		target : node
			The target node.
		storePred : bool
			If true, the algorithm will also store the predecessors
			and reconstruct a shortest path from @a source and @a target.
	"""
		if type(self) == STSP:
			raise RuntimeError("Error, you may not use STSP directly, use a sub-class instead")

	def getPath(self):
		"""
		Returns a shortest path from the source node to the target node (without
		including them). Note: the shortest path can be constructed only if the
		algorithm is executed with @a storePred set to true.

		Returns:
		--------
		vector
			A shortest path from the source node to the target node.
		"""
		return (<_STSP*>(self._this)).getPath()

	def getPredecessors(self):
		"""
		Returns the predecessor nodes from the target node to the source node,
		Note: predecessors are stored only if the algorithm is executed with
		storePred set to true.

		Returns:
		--------
		vector
			The list of predecessors from @a target to @a source.
		"""
		return (<_STSP*>(self._this)).getPredecessors()

	def getDistance(self):
		"""
		Returns the distance from the source node to the target node

		Returns:
		--------
		edgeweight
			The distance from source to the target node.
		"""
		return (<_STSP*>(self._this)).getDistance()

cdef extern from "<networkit/distance/SSSP.hpp>":

	cdef cppclass _SSSP "NetworKit::SSSP"(_Algorithm):
		_SSSP(_Graph G, node source, bool_t storePaths, bool_t storeNodesSortedByDistance, node target) except +
		vector[edgeweight] getDistances() except +
		edgeweight distance(node t) except +
		vector[node] getPredecessors(node t) except +
		vector[node] getPath(node t, bool_t forward) except +
		set[vector[node]] getPaths(node t, bool_t forward) except +
		vector[node] getNodesSortedByDistance() except +
		double _numberOfPaths(node t) except +
		void setSource(node s) except +
		void setTarget(node t) except +

cdef class SSSP(Algorithm):
	""" Base class for single source shortest path algorithms. """
	cdef Graph _G

	def __init__(self, *args, **namedargs):
		if type(self) == SSSP:
			raise RuntimeError("Error, you may not use SSSP directly, use a sub-class instead")

	def getDistances(self):
		"""
		Returns a list of weighted distances from the source node, i.e. the
 	 	length of the shortest path from the source node to any other node.

 	 	Returns:
 	 	--------
 	 	vector
 	 		The weighted distances from the source node to any other node in the graph.
		"""
		return (<_SSSP*>(self._this)).getDistances()

	def distance(self, t):
		"""
		Returns the distance from the source node to @a t.

		Parameters:
		-----------
		t : node
			Target node.

		Returns:
		--------
		double
			Distance from the source node to @a t.
		"""
		return (<_SSSP*>(self._this)).distance(t)

	def getPredecessors(self, t):
		"""
		Returns the predecessor nodes of @a t on all shortest paths from source
		to @a t.
		Parameters:
		-----------
		t : node
			Target node.

		Returns:
		--------
		list
			The predecessors of @a t on all shortest paths from source to @a t.
		"""
		return (<_SSSP*>(self._this)).getPredecessors(t)

	def getPath(self, t, forward=True):
		"""
		Returns a shortest path from source to @a t and an empty path if source and @a t
		are not connected.

		Parameters:
		-----------
		t : node
			Target node.
		forward : bool
			If @c true (default) the path is directed from source to @a t, otherwise the path
			is reversed.

		Returns:
		--------
		list
			A shortest path from source to @a t or an empty path.
		"""
		return (<_SSSP*>(self._this)).getPath(t, forward)

	def getPaths(self, t, forward=True):
		"""
		Returns all shortest paths from source to @a t and an empty set if source
		and @a t are not connected.

		Parameters:
		-----------
		t : node
			Target node.
		forward : bool
			If @c true (default) the path is directed from source to
			@a t, otherwise the path is reversed.

		Returns:
		--------
			All shortest paths from source node to target node @a t.
		"""
		cdef set[vector[node]] paths = (<_SSSP*>(self._this)).getPaths(t, forward)
		result = []
		for elem in paths:
			result.append(list(elem))
		return result

	def getNodesSortedByDistance(self):
		""" Returns a list of nodes ordered in increasing distance from the source.

		For this functionality to be available, storeNodesSortedByDistance has to be set to true in the constructor.
		There are no guarantees regarding the ordering of two nodes with the same distance to the source.

		Returns:
		--------
		list
			Nodes ordered in increasing distance from the source.
		"""
		return (<_SSSP*>(self._this)).getNodesSortedByDistance()

	def numberOfPaths(self, t):
		"""
		Returns the number of paths from the source node to @a t.

		Parameters:
		-----------
		t : node
			Target node.

		Returns:
		--------
		int
			The number of paths from the source node to @a t.
		"""
		return (<_SSSP*>(self._this))._numberOfPaths(t)

	def setSource(self, s not None):
		"""
		Sets a new source node.

		Parameters:
		-----------
		s : node
			New source node.
		"""
		(<_SSSP*>(self._this)).setSource(s)

	def setTarget(self, t not None):
		"""
		Sets a new target node.

		Parameters:
		-----------
		t : node
			New target node.
		"""
		(<_SSSP*>(self._this)).setTarget(t)

cdef extern from "<networkit/distance/DynSSSP.hpp>":

	cdef cppclass _DynSSSP "NetworKit::DynSSSP"(_SSSP):
		_DynSSSP(_Graph G, node source, bool_t storePaths, bool_t storeStack, node target) except +
		void update(_GraphEvent ev) except +
		void updateBatch(vector[_GraphEvent] batch) except +
		bool_t modified() except +
		void setTargetNode(node t) except +

cdef class DynSSSP(SSSP):
	""" Base class for single source shortest path algorithms in dynamic graphs. """
	def __init__(self, *args, **namedargs):
		if type(self) == SSSP:
			raise RuntimeError("Error, you may not use DynSSSP directly, use a sub-class instead")

	def update(self, ev):
		""" Updates shortest paths with the edge insertion.

		Parameters:
		-----------
		ev : GraphEvent.
		"""
		(<_DynSSSP*>(self._this)).update(_GraphEvent(ev.type, ev.u, ev.v, ev.w))

	def updateBatch(self, batch):
		""" Updates shortest paths with the batch `batch` of edge insertions.

		Parameters:
		-----------
		batch : list of GraphEvent.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		(<_DynSSSP*>(self._this)).updateBatch(_batch)

	def modified(self):
		return (<_DynSSSP*>(self._this)).modified()

	def setTargetNode(self, t):
		(<_DynSSSP*>(self._this)).setTargetNode(t)

cdef extern from "<networkit/distance/AdamicAdarDistance.hpp>":

	cdef cppclass _AdamicAdarDistance "NetworKit::AdamicAdarDistance":
		_AdamicAdarDistance(const _Graph& G) except +
		void preprocess() except +
		double distance(node u, node v) except +
		vector[double] getEdgeScores() except +

cdef class AdamicAdarDistance:
	"""
	Calculate the adamic adar similarity.

	Parameters:
	-----------
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
		Returns:
		--------
		vector[double]
			The edge attribute that contains the adamic adar similarity.

		"""
		#### TODO: convert distance to similarity!?! ####
		return self._this.getEdgeScores()

	def distance(self, node u, node v):
		return self._this.distance(u, v)

cdef extern from "<networkit/distance/Diameter.hpp>" namespace "NetworKit":

	cdef enum _DiameterAlgo "NetworKit::DiameterAlgo":
		automatic = 0
		exact = 1
		estimatedRange = 2
		estimatedSamples = 3
		estimatedPedantic = 4

class DiameterAlgo(object):
	Automatic = automatic
	Exact = exact
	EstimatedRange = estimatedRange
	EstimatedSamples = estimatedSamples
	EstimatedPedantic = estimatedPedantic

cdef extern from "<networkit/distance/Diameter.hpp>" namespace "NetworKit::Diameter":

	cdef cppclass _Diameter "NetworKit::Diameter"(_Algorithm):
		_Diameter(_Graph G, _DiameterAlgo algo, double error, count nSamples) except +
		pair[count, count] getDiameter() nogil except +

cdef class Diameter(Algorithm):
	cdef Graph _G
	"""
	TODO: docstring
	"""
	def __cinit__(self, Graph G not None, algo = DiameterAlgo.Automatic, error = -1., nSamples = 0):
		self._G = G
		self._this = new _Diameter(G._this, algo, error, nSamples)

	def getDiameter(self):
		return (<_Diameter*>(self._this)).getDiameter()

cdef extern from "<networkit/distance/Eccentricity.hpp>" namespace "NetworKit::Eccentricity":

	pair[node, count] getValue(_Graph G, node v) except +

cdef class Eccentricity:
	"""
	The eccentricity of a node `u` is defined as the distance to the farthest node from node u. In other words, it is the longest shortest-path starting from node `u`.
	"""

	@staticmethod
	def getValue(Graph G, v):
		"""
		Returns:
		--------
		pair[node, count]
			node is the farthest node `v` from `u`, and the count is the length of the shortest path from `u` to `v`.
		"""
		return getValue(G._this, v)

cdef extern from "<networkit/distance/EffectiveDiameterApproximation.hpp>" namespace "NetworKit::EffectiveDiameterApproximation":

	cdef cppclass _EffectiveDiameterApproximation "NetworKit::EffectiveDiameterApproximation"(_Algorithm):
		_EffectiveDiameterApproximation(_Graph& G, double ratio, count k, count r) except +
		double getEffectiveDiameter() except +

cdef class EffectiveDiameterApproximation(Algorithm):
	"""
	Calculates the effective diameter of a graph.
	The effective diameter is defined as the number of edges on average to reach a given ratio of all other nodes.

	Implementation after the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining in Massive Graphs"[1]

	[1] by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	ratio : double
		The percentage of nodes that shall be within stepwidth, default = 0.9
	k : count
		number of parallel approximations, bigger k -> longer runtime, more precise result; default = 64
	r : count
		number of additional bits, important in tiny graphs; default = 7
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, double ratio=0.9, count k=64, count r=7):
		self._G = G
		self._this = new _EffectiveDiameterApproximation(G._this, ratio, k, r)

	def getEffectiveDiameter(self):
		"""
		Returns:
		--------
		double
			the approximated effective diameter
		"""
		return (<_EffectiveDiameterApproximation*>(self._this)).getEffectiveDiameter()

cdef extern from "<networkit/distance/EffectiveDiameter.hpp>" namespace "NetworKit::EffectiveDiameter":

	cdef cppclass _EffectiveDiameter "NetworKit::EffectiveDiameter"(_Algorithm):
		_EffectiveDiameter(_Graph& G, double ratio) except +
		double getEffectiveDiameter() except +

cdef class EffectiveDiameter(Algorithm):
	"""
	Calculates the effective diameter of a graph.
	The effective diameter is defined as the number of edges on average to reach a given ratio of all other nodes.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	ratio : double
		The percentage of nodes that shall be within stepwidth; default = 0.9
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, double ratio=0.9):
		self._G = G
		self._this = new _EffectiveDiameter(G._this, ratio)

	def getEffectiveDiameter(self):
		"""
		Returns:
		--------
		double
			the effective diameter
		"""
		return (<_EffectiveDiameter*>(self._this)).getEffectiveDiameter()

cdef extern from "<networkit/distance/HopPlotApproximation.hpp>" namespace "NetworKit::HopPlotApproximation":

	cdef cppclass _HopPlotApproximation "NetworKit::HopPlotApproximation"(_Algorithm):
		_HopPlotApproximation(_Graph& G, count maxDistance, count k, count r) except +
		map[count, double] getHopPlot() except +

cdef class HopPlotApproximation(Algorithm):
	"""
	Computes an approxmation of the hop-plot of a given graph.
	The hop-plot is the set of pairs (d, g(g)) for each natural number d
	and where g(d) is the fraction of connected node pairs whose shortest connecting path has length at most d.

	Implementation after the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining in Massive Graphs"[1]

	[1] by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	maxDistance : double
		maximum distance between considered nodes
		set to 0 or negative to get the hop-plot for the entire graph so that each node can reach each other node
	k : count
		number of parallel approximations, bigger k -> longer runtime, more precise result; default = 64
	r : count
		number of additional bits, important in tiny graphs; default = 7
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, count maxDistance=0, count k=64, count r=7):
		self._G = G
		self._this = new _HopPlotApproximation(G._this, maxDistance, k, r)

	def getHopPlot(self):
		"""
		Returns:
		--------
		map
			number of connected nodes for each distance
		"""
		cdef map[count, double] hp = (<_HopPlotApproximation*>(self._this)).getHopPlot()
		result = dict()
		for elem in hp:
			result[elem.first] = elem.second
		return result

cdef extern from "<networkit/distance/NeighborhoodFunction.hpp>" namespace "NetworKit::NeighborhoodFunction":

	cdef cppclass _NeighborhoodFunction "NetworKit::NeighborhoodFunction"(_Algorithm):
		_NeighborhoodFunction(_Graph& G) except +
		vector[count] getNeighborhoodFunction() except +

cdef class NeighborhoodFunction(Algorithm):
	"""
	Computes the neighborhood function exactly.
	The neighborhood function N of a graph G for a given distance t is defined
	as the number of node pairs (u,v) that can be reached within distance t.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None):
		self._G = G
		self._this = new _NeighborhoodFunction(G._this)

	def getNeighborhoodFunction(self):
		"""
		Returns:
		--------
		list
			the i-th element denotes the number of node pairs that have a distance at most (i+1)
		"""
		return (<_NeighborhoodFunction*>(self._this)).getNeighborhoodFunction()

cdef extern from "<networkit/distance/NeighborhoodFunctionApproximation.hpp>" namespace "NetworKit::NeighborhoodFunctionApproximation":

	cdef cppclass _NeighborhoodFunctionApproximation "NetworKit::NeighborhoodFunctionApproximation"(_Algorithm):
		_NeighborhoodFunctionApproximation(_Graph& G, count k, count r) except +
		vector[count] getNeighborhoodFunction() except +

cdef class NeighborhoodFunctionApproximation(Algorithm):
	"""
	Computes an approximation of the neighborhood function.
	The neighborhood function N of a graph G for a given distance t is defined
	as the number of node pairs (u,v) that can be reached within distance t.

	Implementation after the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining in Massive Graphs"[1]

	[1] by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	k : count
		number of approximations, bigger k -> longer runtime, more precise result; default = 64
	r : count
		number of additional bits, important in tiny graphs; default = 7
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, count k=64, count r=7):
		self._G = G
		self._this = new _NeighborhoodFunctionApproximation(G._this, k, r)

	def getNeighborhoodFunction(self):
		"""
		Returns:
		--------
		list
			the i-th element denotes the number of node pairs that have a distance at most (i+1)
		"""
		return (<_NeighborhoodFunctionApproximation*>(self._this)).getNeighborhoodFunction()

cdef extern from "<networkit/distance/Volume.hpp>" namespace "NetworKit::Volume":

	double volume(const _Graph G, const double r, const count samples) nogil except +
	vector[double] volume(const _Graph G, const vector[double] r, const count samples) nogil except +

cdef class Volume:

	@staticmethod
	def volume(Graph G, r, count samples=500):
		"""
		Number of nodes within a given radius (or radii); average for many nodes

		Please find further information about the volume and its meaning in the
		following publication:

		Franz-Benjamin Mocnik: "The Polynomial Volume Law of Complex Networks in
		the Context of Local and Global Optimization", Scientific Reports 8(11274)
		2018. doi: 10.1038/s41598-018-29131-0

		Parameters:
		-----------
		G : networkit.Graph
			the graph
		r : double
			the radius (or radii)
		samples : count
			the number of samples
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
		_JaccardDistance(const _Graph& G, const vector[count]& triangles) except +
		void preprocess() except +
		vector[double] getEdgeScores() except +

cdef class JaccardDistance:
	"""
	The Jaccard distance measure assigns to each edge the jaccard coefficient
	of the neighborhoods of the two adjacent nodes.

	Parameters:
	-----------
	G : networkit.Graph
		The graph to calculate Jaccard distances for.
	triangles : vector[count]
		Previously calculated edge triangle counts.
	"""

	cdef _JaccardDistance* _this
	cdef Graph _G
	cdef vector[count] triangles

	def __cinit__(self, Graph G, vector[count] triangles):
		self._G = G
		self._triangles = triangles
		self._this = new _JaccardDistance(G._this, self._triangles)

	def __dealloc__(self):
		del self._this

	def getAttribute(self):
		return self._this.getEdgeScores()

cdef class JaccardSimilarityAttributizer:
	"""
	The Jaccard similarity measure assigns to each edge (1 - the jaccard coefficient
	of the neighborhoods of the two adjacent nodes).

	Parameters:
	-----------
	G : networkit.Graph
		The graph to calculate Jaccard similarities for.
	triangles : vector[count]
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
		#convert distance to similarity
		self._this.preprocess()
		return [1 - x for x in self._this.getEdgeScores()]

cdef extern from "<networkit/distance/AlgebraicDistance.hpp>":

	cdef cppclass _AlgebraicDistance "NetworKit::AlgebraicDistance":
		_AlgebraicDistance(_Graph G, count numberSystems, count numberIterations, double omega, index norm, bool_t withEdgeScores) except +
		void preprocess() except +
		double distance(node, node) except +
		vector[double] getEdgeScores() except +


cdef class AlgebraicDistance:
	"""
	Algebraic distance assigns a distance value to pairs of nodes
    according to their structural closeness in the graph.
    Algebraic distances will become small within dense subgraphs.

	Parameters:
	-----------
	G : networkit.Graph
		The graph to calculate Jaccard distances for.
	numberSystems : count
	 	Number of vectors/systems used for algebraic iteration.
	numberIterations : count
	 	Number of iterations in each system.
	omega : double
	 	attenuation factor in [0,1] influencing convergence speed.
	norm : index
		The norm factor of the extended algebraic distance.
	withEdgeScores : bool
		calculate array of scores for edges {u,v} that equal ad(u,v)
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
	""" Computes the Euclidean Commute Time Distance between each pair of nodes for an undirected unweighted graph.

	CommuteTimeDistance(G)

	Create CommuteTimeDistance for Graph `G`.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	tol: double
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G, double tol = 0.1):
		self._G = G
		self._this = new _CommuteTimeDistance(G._this, tol)

	def runApproximation(self):
		""" Computes approximation of the ECTD. """
		return (<_CommuteTimeDistance*>(self._this)).runApproximation()

	def runParallelApproximation(self):
		""" Computes approximation (in parallel) of the ECTD. """
		return (<_CommuteTimeDistance*>(self._this)).runParallelApproximation()

	def distance(self, u, v):
		"""  Returns the ECTD between node u and node v.

		u : node
		v : node
		"""
		return (<_CommuteTimeDistance*>(self._this)).distance(u, v)

	def runSinglePair(self, u, v):
		"""  Returns the ECTD between node u and node v, without preprocessing.

		u : node
		v : node
		"""
		return (<_CommuteTimeDistance*>(self._this)).runSinglePair(u, v)

	def runSingleSource(self, u):
		"""  Returns the sum of the ECTDs from u, without preprocessing.

		u : node
		"""
		return (<_CommuteTimeDistance*>(self._this)).runSingleSource(u)

cdef extern from "<networkit/distance/NeighborhoodFunctionHeuristic.hpp>" namespace "NetworKit::NeighborhoodFunctionHeuristic::SelectionStrategy":

	enum _SelectionStrategy "NetworKit::NeighborhoodFunctionHeuristic::SelectionStrategy":
		RANDOM
		SPLIT

cdef extern from "<networkit/distance/NeighborhoodFunctionHeuristic.hpp>" namespace "NetworKit::NeighborhoodFunctionHeuristic":

	cdef cppclass _NeighborhoodFunctionHeuristic "NetworKit::NeighborhoodFunctionHeuristic"(_Algorithm):
		_NeighborhoodFunctionHeuristic(_Graph& G, const count nSamples, const _SelectionStrategy strategy) except +
		vector[count] getNeighborhoodFunction() except +

cdef class NeighborhoodFunctionHeuristic(Algorithm):
	"""
	Computes a heuristic of the neighborhood function.
	The algorithm runs nSamples breadth-first searches and scales the results up to the actual amount of nodes.
	Accepted strategies are "split" and "random".

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	nSamples : count
		the amount of samples, set to zero for heuristic of max(sqrt(m), 0.15*n)
	strategy : enum
		the strategy to select the samples, accepts "random" or "split"
	"""
	cdef Graph _G

	RANDOM = 0
	SPLIT = 1

	def __cinit__(self, Graph G not None, count nSamples=0, strategy=SPLIT):
		self._G = G
		self._this = new _NeighborhoodFunctionHeuristic(G._this, nSamples, strategy)

	def getNeighborhoodFunction(self):
		"""
		Returns:
		--------
		list
			the i-th element denotes the number of node pairs that have a distance at most (i+1)
		"""
		return (<_NeighborhoodFunctionHeuristic*>(self._this)).getNeighborhoodFunction()

cdef extern from "<networkit/distance/APSP.hpp>":

	cdef cppclass _APSP "NetworKit::APSP"(_Algorithm):
		_APSP(_Graph G) except +
		vector[vector[edgeweight]] getDistances() except +
		edgeweight getDistance(node u, node v) except +

cdef class APSP(Algorithm):
	""" All-Pairs Shortest-Paths algorithm (implemented running Dijkstra's algorithm from each node, or BFS if G is unweighted).

    APSP(G)

    Computes all pairwise shortest-path distances in G.

    Parameters:
	-----------
	G : networkit.Graph
		The graph.
    """
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _APSP(G._this)

	def __dealloc__(self):
		self._G = None

	def getDistances(self):
		""" Returns a vector of vectors of distances between each node pair.

 	 	Returns:
 	 	--------
 	 	vector of vectors
 	 		The shortest-path distances from each node to any other node in the graph.
		"""
		return (<_APSP*>(self._this)).getDistances()

	def getDistance(self, node u, node v):
		""" Returns the length of the shortest path from source 'u' to target 'v'.

		Parameters:
		-----------
		u : node
			Source node.
		v : node
			Target node.

		Returns:
		--------
		int or float
			The distance from 'u' to 'v'.
		"""
		return (<_APSP*>(self._this)).getDistance(u, v)

cdef extern from "<networkit/distance/SPSP.hpp>":

	cdef cppclass _SPSP "NetworKit::SPSP"(_Algorithm):
		_SPSP(_Graph G, vector[node].iterator sourcesFirst, vector[node].iterator sourcesLast) except +
		vector[vector[edgeweight]] getDistances() except +
		edgeweight getDistance(node u, node v) except +
		void setSources(vector[node].iterator sourcesFirst, vector[node].iterator sourcesLast)

cdef class SPSP(Algorithm):
	""" Some-Pairs Shortest-Paths algorithm (implemented running Dijkstra's algorithm from each source
		node, or BFS if G is unweighted).

    SPSP(G)

    Computes pairwise shortest-path distances from the source nodes to all the nodes in G.

    Parameters:
	-----------
	G : networkit.Graph
		The graph.
	sources : list.
		Set of source nodes.
    """
	cdef Graph _G

	def __cinit__(self, Graph G not None, vector[node] sources):
		self._G = G
		self._this = new _SPSP(G._this, sources.begin(), sources.end())

	def __dealloc__(self):
		self._G = None

	def getDistances(self):
		""" Returns a vector of vectors of distances between each source node
		and all the other nodes pair.

 	 	Returns:
 	 	--------
 	 	vector of vectors
			The shortest-path distances from each source node to any other node
			in the graph.
		"""
		return (<_SPSP*>self._this).getDistances()

	def getDistance(self, node u, node v):
		""" Returns the length of the shortest path from source 'u' to target 'v'.

		Parameters:
		-----------
		u : node
			Source node.
		v : node
			Target node.

		Returns:
		--------
		float
			The distance from 'u' to 'v'.
		"""
		return (<_SPSP*>self._this).getDistance(u, v)

	def setSources(self, vector[node] sources):
		""" Sets the source nodes.

		Parameters
		----------
		sources : list
			List of the new source nodes.
		"""
		(<_SPSP*>self._this).setSources(sources.begin(), sources.end())


cdef extern from "<networkit/distance/DynAPSP.hpp>":

	cdef cppclass _DynAPSP "NetworKit::DynAPSP"(_APSP):
		_DynAPSP(_Graph G) except +
		void update(_GraphEvent ev) except +
		void updateBatch(vector[_GraphEvent] batch) except +

cdef class DynAPSP(APSP):
	""" All-Pairs Shortest-Paths algorithm for dynamic graphs.

	DynAPSP(G)

	Computes all pairwise shortest-path distances in G.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
		"""
	def __init__(self, Graph G):
		self._G = G
		self._this = new _DynAPSP(G._this)

	def update(self, ev):
		""" Updates shortest paths with the edge insertion.

		Parameters:
		-----------
		ev : GraphEvent.
		"""
		(<_DynAPSP*>(self._this)).update(_GraphEvent(ev.type, ev.u, ev.v, ev.w))

	def updateBatch(self, batch):
		""" Updates shortest paths with the batch `batch` of edge insertions.

		Parameters:
		-----------
		batch : list of GraphEvent.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		(<_DynAPSP*>(self._this)).updateBatch(_batch)

cdef extern from "<networkit/distance/BFS.hpp>":

	cdef cppclass _BFS "NetworKit::BFS"(_SSSP):
		_BFS(_Graph G, node source, bool_t storePaths, bool_t storeNodesSortedByDistance, node target) except +

cdef class BFS(SSSP):
	""" Simple breadth-first search on a Graph from a given source

	BFS(G, source, storePaths=True, storeNodesSortedByDistance=False, target=None)

	Create BFS for `G` and source node `source`.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	source : node
		The source node of the breadth-first search.
	storePaths : bool
		store paths and number of paths?
	target: node
		terminate search when the target has been reached
	"""

	def __cinit__(self, Graph G, source, storePaths=True, storeNodesSortedByDistance=False, target=none):
		self._G = G
		self._this = new _BFS(G._this, source, storePaths, storeNodesSortedByDistance, target)


cdef extern from "<networkit/distance/Dijkstra.hpp>":

	cdef cppclass _Dijkstra "NetworKit::Dijkstra"(_SSSP):
		_Dijkstra(_Graph G, node source, bool_t storePaths, bool_t storeNodesSortedByDistance, node target) except +

cdef class Dijkstra(SSSP):
	""" Dijkstra's SSSP algorithm.
	Returns list of weighted distances from node source, i.e. the length of the shortest path from source to
	any other node.

    Dijkstra(G, source, storePaths=True, storeNodesSortedByDistance=False, target=None)

    Creates Dijkstra for `G` and source node `source`.

    Parameters:
	-----------
	G : networkit.Graph
		The graph.
	source : node
		The source node.
	storePaths : bool
		Paths are reconstructable and the number of paths is stored.
	storeNodesSortedByDistance: bool
		Store a vector of nodes ordered in increasing distance from the source.
	target : node
		target node. Search ends when target node is reached. t is set to None by default.
    """
	def __cinit__(self, Graph G, source, storePaths=True, storeNodesSortedByDistance=False, node target=none):
		self._G = G
		self._this = new _Dijkstra(G._this, source, storePaths, storeNodesSortedByDistance, target)

cdef extern from "<networkit/distance/DynBFS.hpp>":

	cdef cppclass _DynBFS "NetworKit::DynBFS"(_DynSSSP):
		_DynBFS(_Graph G, node source) except +

cdef class DynBFS(DynSSSP):
	""" Dynamic version of BFS.

	DynBFS(G, source)

	Create DynBFS for `G` and source node `source`.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	source : node
		The source node of the breadth-first search.
	storeStack : bool
		maintain a stack of nodes in order of decreasing distance?
	"""
	def __cinit__(self, Graph G, source):
		self._G = G
		self._this = new _DynBFS(G._this, source)

cdef extern from "<networkit/distance/DynDijkstra.hpp>":

	cdef cppclass _DynDijkstra "NetworKit::DynDijkstra"(_DynSSSP):
		_DynDijkstra(_Graph G, node source) except +

cdef class DynDijkstra(DynSSSP):
	""" Dynamic version of Dijkstra.

	DynDijkstra(G, source)

	Create DynDijkstra for `G` and source node `source`.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	source : node
		The source node of the breadth-first search.

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
		count getHops() except +

cdef class BidirectionalBFS(STSP):
	"""
		Implements a bidirectional breadth-first search on a graph from
		two given source and target nodes.
		Explores the graph from both the source and target nodes until
		the two explorations meet.

		Parameters:
		-----------

		G : networkit.Graph
			The input graph.
		source : node
			The source node.
		target : node
			The target node.
		storePred : bool
			If true, the algorithm will also store the predecessors
			and reconstruct a shortest path from @a source and @a target.
	"""

	def __cinit__(self, Graph G, node source, node target, bool_t storePred=True):
		self._this = new _BidirectionalBFS(G._this, source, target, storePred)

	def getHops(self):
		"""
		Returns the distance (i.e., number of hops) from the source to the
		target node.

		Returns:
		--------
		count
			Number of hops from the source to the target node.
		"""
		return (<_BidirectionalBFS*>(self._this)).getHops()

cdef extern from "<networkit/distance/BidirectionalDijkstra.hpp>":
	cdef cppclass _BidirectionalDijkstra "NetworKit::BidirectionalDijkstra"(_STSP):
		_BidirectionalDijkstra(_Graph G, node source, node target, bool_t storePred) except +

cdef class BidirectionalDijkstra(STSP):
	"""
		Bidirectional implementation of the Dijkstra algorithm from
		two given source and target nodes.
		Explores the graph from both the source and target nodes until
		the two explorations meet.

		Parameters:
		-----------

		G : networkit.Graph
			The input graph.
		source : node
			The source node.
		target : node
			The target node.
		storePred : bool
			If true, the algorithm will also store the predecessors
			and reconstruct a shortest path from @a source and @a target.
	"""

	def __cinit__(self, Graph G, node source, node target, bool_t storePred=True):
		self._this = new _BidirectionalDijkstra(G._this, source, target, storePred)

cdef extern from "<networkit/distance/AStar.hpp>":
	cdef cppclass _AStar "NetworKit::AStar"(_STSP):
		_AStar(_Graph G, vector[double] &heu, node source, node target, bool_t storePred) except +

cdef class AStar(STSP):
	"""
	A* path-finding algorithm.

	Parameters:
	-----------

	G : networkit.Graph
		The input graph.
	heu : list
		List of lower bounds of the distance of each node to the target.
	source : node
		The source node.
	target : node
		The target node.
	storePred : bool
		If true, the algorithm will also store the predecessors
		and reconstruct a shortest path from @a source and @a target.
	"""

	cdef vector[double] heu
	def __cinit__(self, Graph G, vector[double] &heu, node source, node target, bool_t storePred=True):
		self.heu = heu
		self._this = new _AStar(G._this, self.heu, source, target, storePred)

cdef extern from "<networkit/distance/AllSimplePaths.hpp>":

	cdef cppclass _AllSimplePaths "NetworKit::AllSimplePaths":
		_AllSimplePaths(_Graph G, node source, node target, count cutoff) except +
		void run() nogil except +
		count numberOfSimplePaths() except +
		vector[vector[node]] getAllSimplePaths() except +
		void forAllSimplePaths[Callback](Callback c) except +

cdef class AllSimplePaths:
	""" Algorithm to compute all existing simple paths from a source node to a target node. The maximum length of the paths can be fixed through 'cutoff'.
		CAUTION: This algorithm could take a lot of time on large networks (many edges), especially if the cutoff value is high or not specified.

	AllSimplePaths(G, source, target, cutoff=none)

	Create AllSimplePaths for `G`, source node `source`, target node 'target' and cutoff 'cutoff'.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	source : node
		The source node.
	target : node
		The target node.
	cutoff : count
		(optional) The maximum length of the simple paths.

	"""

	cdef _AllSimplePaths* _this
	cdef Graph _G

	def __cinit__(self,  Graph G, source, target, cutoff=none):
		self._G = G
		self._this = new _AllSimplePaths(G._this, source, target, cutoff)

	def __dealloc__(self):
		del self._this

	def run(self):
		self._this.run()
		return self

	def numberOfSimplePaths(self):
		"""
		Returns the number of simple paths.

		Returns:
		--------
		count
			The number of simple paths.
		"""
		return self._this.numberOfSimplePaths()

	def getAllSimplePaths(self):
		"""
		Returns all the simple paths from source to target.

		Returns:
		--------
		A vector of vectors.
			A vector containing vectors which represent all simple paths.
		"""
		return self._this.getAllSimplePaths()

	def forAllSimplePaths(self, object callback):
		""" More efficient path iterator. Iterates over all the simple paths.

		Parameters:
		-----------
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
	""" Simple reverse breadth-first search on a Graph from a given source

	ReverseBFS(G, source, storePaths=True, storeNodesSortedByDistance=False, target=None)

	Create ReverseBFS for `G` and source node `source`.

	Parameters:
	-----------
	G : networkit.Graph
		The graph.
	source : node
		The source node of the breadth-first search.
	storePaths : bool
		Paths are reconstructable and the number of paths is stored.
	target: node
		terminate search when the target has been reached
	"""

	def __cinit__(self, Graph G, source, storePaths=True, storeNodesSortedByDistance=False, target=none):
		self._G = G
		self._this = new _ReverseBFS(G._this, source, storePaths, storeNodesSortedByDistance, target)

