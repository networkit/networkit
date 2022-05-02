# distutils: language=c++

from libc.stdint cimport uint64_t
from libc.stdint cimport uint8_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp cimport bool as bool_t

import math
from enum import Enum

from .base cimport _Algorithm, Algorithm
from .dynamics cimport _GraphEvent, GraphEvent
from .graph cimport _Graph, Graph
from .structures cimport _Cover, Cover, _Partition, Partition, count, index, node, edgeweight
from networkit.algebraic import adjacencyEigenvector, PageRankMatrix, symmetricEigenvectors

cdef extern from "limits.h":
	cdef uint64_t ULONG_MAX

cdef extern from "<networkit/centrality/Centrality.hpp>":

	cdef cppclass _Centrality "NetworKit::Centrality"(_Algorithm):
		_Centrality(_Graph, bool_t, bool_t) except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +
		double maximum() except +
		double centralization() except +

cdef class Centrality(Algorithm):
	""" Abstract base class for centrality measures"""

	cdef Graph _G

	def __init__(self, *args, **kwargs):
		if type(self) == Centrality:
			raise RuntimeError("Error, you may not use Centrality directly, use a sub-class instead")

	def __dealloc__(self):
		self._G = None # just to be sure the graph is deleted

	def scores(self):
		"""
		scores()

		Returns the scores of all nodes for the centrality algorithm.

		Returns
		-------
		list(float)
			The list of all scores.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).scores()

	def score(self, v):
		"""
		score(v)
		
		Returns the score of node v for the centrality algorithm.

		Parameters
		----------
		v : int
			Node index.

		Returns
		-------
		float
			The score of node v.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).score(v)

	def ranking(self):
		"""
		ranking()

		Returns the ranking of nodes according to their score.

		Returns
		-------
		dict(tuple(int, float))
			A vector of pairs sorted into descending order. Each pair contains a node and the corresponding score
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).ranking()

	def maximum(self):
		"""
		maximum()

		Return the maximum theoretical centrality score.

		Returns
		-------
		float
			The maximum theoretical centrality score for the given graph.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).maximum()

	def centralization(self):
		"""
		centralization()

		Compute the centralization of a network with respect to some centrality measure.
		The centralization of any network is a measure of how central its most central
		node is in relation to how central all the other nodes are.
		Centralization measures then (a) calculate the sum in differences
		in centrality between the most central node in a network and all other nodes;
		and (b) divide this quantity by the theoretically largest such sum of
		differences in any network of the same size.

		Returns
		-------
		float
			Centralization value.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_Centrality*>(self._this)).centralization()

cdef extern from "<networkit/centrality/Betweenness.hpp>":

	cdef cppclass _Betweenness "NetworKit::Betweenness" (_Centrality):
		_Betweenness(_Graph, bool_t, bool_t) except +
		vector[double] edgeScores() except +

cdef class Betweenness(Centrality):
	"""
	Betweenness(G, normalized=False, computeEdgeCentrality=False)

	Constructs the Betweenness class for the given Graph `G`. If the betweenness scores should be normalized,
	then set `normalized` to True. The run() method takes O(nm) time, where n is the number
 	of nodes and m is the number of edges of the graph.

 	Parameters
 	----------
 	G : networkit.Graph
 		The input graph.
 	normalized : bool, optional
 		Set this parameter to True if scores should be normalized in the interval [0,1]. Default: False
	computeEdgeCentrality: bool, optional
		Set this to true if edge betweenness scores should be computed as well. Default: False
	"""

	def __cinit__(self, Graph G, normalized=False, computeEdgeCentrality=False):
		self._G = G
		self._this = new _Betweenness(G._this, normalized, computeEdgeCentrality)


	def edgeScores(self):
		"""
		edgeScores()
		
		Get a vector containing the betweenness score for each edge in the graph
		in ascending edge ID order.

		Returns
		-------
		list(float)
			The betweenness scores calculated by run().
		"""
		return (<_Betweenness*>(self._this)).edgeScores()

cdef extern from "<networkit/centrality/ApproxBetweenness.hpp>":

	cdef cppclass _ApproxBetweenness "NetworKit::ApproxBetweenness" (_Centrality):
		_ApproxBetweenness(_Graph, double, double, double) except +
		count numberOfSamples() except +

cdef class ApproxBetweenness(Centrality):
	"""
 	ApproxBetweenness(G, epsilon=0.01, delta=0.1, universalConstant=1.0)

	Approximation of betweenness centrality according to algorithm described in Matteo Riondato 
	and Evgenios M. Kornaropoulos: Fast Approximation of Betweenness Centrality through Sampling

 	The algorithm approximates the betweenness of all vertices so that the scores are
	within an additive error epsilon with probability at least :math:`(1- delta)`.
	The values are normalized by default. The run() method takes O(m) time per sample, where  m is
	the number of edges of the graph. The number of samples is proportional to universalConstant/epsilon^2.
	Although this algorithm has a theoretical guarantee, the algorithm implemented in Estimate Betweenness usually performs better in practice.
	Therefore, we recommend to use EstimateBetweenness if no theoretical guarantee is needed.

	Parameters
	----------
	G : networkit.Graph
		The input graph
	epsilon : double, optional
		Maximum additive error
	delta : double, optional
		Probability that the values are within the error guarantee
	universalConstant: double, optional
		The universal constant to be used in computing the sample size.
		It is 1 by default. Some references suggest using 0.5, but there
		is no guarantee in this case.
	"""

	def __cinit__(self, Graph G, epsilon=0.01, delta=0.1, universalConstant=1.0):
		self._G = G
		self._this = new _ApproxBetweenness(G._this, epsilon, delta, universalConstant)

	def numberOfSamples(self):
		return (<_ApproxBetweenness*>(self._this)).numberOfSamples()


cdef extern from "<networkit/centrality/EstimateBetweenness.hpp>":

	cdef cppclass _EstimateBetweenness"NetworKit::EstimateBetweenness" (_Centrality):
		_EstimateBetweenness(_Graph, count, bool_t, bool_t) except +

cdef class EstimateBetweenness(Centrality):
	""" 
	EstimateBetweenness(G, nSamples, normalized=False, parallel=False)
	
	Estimation of betweenness centrality according to algorithm described in
	Sanders, Geisberger, Schultes: Better Approximation of Betweenness Centrality

	The algorithm estimates the betweenness of all nodes, using weighting
	of the contributions to avoid biased estimation. The run() method takes O(m)
	time per sample, where  m is the number of edges of the graph. There is no proven
	theoretical guarantee on the quality of the approximation. However, the algorithm
	was shown to perform well in practice. If a guarantee is required, use ApproxBetweenness.

	Parameters
	----------
	G : networkit.Graph
		Input graph
	nSamples : count
		User defined number of samples
	normalized : bool, optional
		Normalize centrality values in interval [0,1]
	parallel : bool, optional
		Run in parallel with additional memory cost z + 3z * t
	"""

	def __cinit__(self, Graph G, nSamples, normalized=False, parallel=False):
		self._G = G
		self._this = new _EstimateBetweenness(G._this, nSamples, normalized, parallel)

cdef extern from "<networkit/centrality/KadabraBetweenness.hpp>":

	cdef cppclass _KadabraBetweenness "NetworKit::KadabraBetweenness" (_Algorithm):
		_KadabraBetweenness(_Graph, double, double, count, count, count) except +
		vector[pair[node, double]] ranking() except +
		vector[node] topkNodesList() except +
		vector[double] topkScoresList() except +
		vector[double] scores() except +
		count getNumberOfIterations() except +
		double getOmega() except +

cdef class KadabraBetweenness(Algorithm):
	"""
	KadabraBetweenness(Graph G, err = 0.01, delta = 0.1, k = 0, unionSample = 0, startFactor = 100

	Approximation of the betweenness centrality and computation of the top-k
	nodes with highest betweenness centrality according to the algorithm
	described in Borassi M. and Natale M. (2016): KADABRA is an ADaptive
	Algorithm for Betweenness via Random Approximation.
	Parallel implementation by Van der Grinten A., Angriman E., and
	Meyerhenke H.: Parallel Adaptive Sampling with almost no
	Synchronization, Euro-Par 2019.
	https://link.springer.com/chapter/10.1007/978-3-030-29400-7_31
	ArXiv pre-print: https://arxiv.org/abs/1903.09422.

	If k = 0 the algorithm approximates the betweenness centrality of all
	vertices of the graph so that the scores are within an additive error
	err with probability at least :math:`(1 - err * delta)`. Otherwise, the algorithm
	computes the exact ranking of the top-k nodes with highest betweenness
	centrality.
	The algorithm relies on an adaptive random sampling technique of shortest
	paths and the number of samples in the worst case is :math:`w = ((log(D - 2) +
	log(2/delta))/err^2` samples, where D is the diameter of the graph.
	Thus, the worst-case performance is :math:`O(w * (|E| + |V|))`, but performs better
	in practice.

	Notes
	-----
	In order to work properly, the Kadabra algorithm requires a random seed
	to be previously set with 'useThreadId' set to True. To do this, call the
	setSeed(<your_seed>, True) fuction within the Random module.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	err : float, optional
		Maximum additive error guaranteed when approximating the
		betweenness centrality of all nodes. Default: 0.01
	delta : float, optional
		Probability that the values of the betweenness centrality are
		within the error guarantee. Default: 0.1
	k : int, optional
		The number of top-k nodes to be computed. Set it to zero to
		approximate the betweenness centrality of all the nodes. Default: 0
	unionSample : int, optional
		Algorithm parameter. Default: 0
	startFactor : int, optional
		Algorithm parameter. Default: 100
	"""

	def __cinit__(self, Graph G, err = 0.01, delta = 0.1, k = 0,
				  unionSample = 0, startFactor = 100):
		self._this = new _KadabraBetweenness(G._this, err, delta, k, unionSample,
										   startFactor)

	def ranking(self):
		"""
		ranking()

		Returns the ranking of the nodes according to their approximated
		betweenness centrality.

		Returns
		-------
		list(int, float)
			A list of pairs (node, betweenness) representing the top-k ranking.
		"""
		return (<_KadabraBetweenness*>(self._this)).ranking()

	def topkNodesList(self):
		"""
		topkNodesList()

		Returns Nodes of the graph sorted by their approximated betweenness
		centrality.

		Returns
		-------
		list(int)
			A list with the top-k nodes with highest approximated betweenness
			centrality.
		"""
		return (<_KadabraBetweenness*>(self._this)).topkNodesList()

	def topkScoresList(self):
		"""
		topkScoresList()

		Returns the sorted list of approximated betweenness centrality scores.

		Returns
		-------
		list(float)
			A list with the top-k scores of the nodes with highest approximated
			betweenness centrality.
		"""
		return (<_KadabraBetweenness*>(self._this)).topkScoresList()

	def scores(self):
		"""
		scores()

		Returns the approximated betweenness centrality score of all the nodes of
		the graph.

		Returns
		-------
		list(float)
			A list with the approximated betweenness centrality score of each node of
			the graph.
		"""
		return (<_KadabraBetweenness*>(self._this)).scores()

	def getNumberOfIterations(self):
		"""
		getNumberOfIterations()

		Returns the total number of samples.

		Returns
		-------
		int
			The total number of shortest paths sampled by the algorithm.
		"""
		return (<_KadabraBetweenness*>(self._this)).getNumberOfIterations()

	def getOmega(self):
		"""
		getOmega()

		Returns the upper bound of the required number of samples.

		Returns
		-------
		int
			Upper bound of the number of shortest paths to be sampled.
		"""
		return(<_KadabraBetweenness*>(self._this)).getOmega()

cdef extern from "<networkit/centrality/DynBetweenness.hpp>":

	cdef cppclass _DynBetweenness "NetworKit::DynBetweenness"(_Algorithm):
		_DynBetweenness(_Graph) except +
		void update(_GraphEvent) except +
		void updateBatch(vector[_GraphEvent]) except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +

cdef class DynBetweenness(Algorithm):
	""" 
	DynBetweenness(G)
	
	The algorithm computes the betweenness centrality of all nodes
	and updates them after an edge insertion.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _DynBetweenness(G._this)

	def update(self, ev):
		"""
		update(ev)
		
		Updates the betweenness centralities after the edge insertions.

		Parameters
		----------
		ev : networkit.dynamics.GraphEvent
			Update the Betweenness based on a given graph event.
		"""
		(<_DynBetweenness*>(self._this)).update(_GraphEvent(ev.type, ev.u, ev.v, ev.w))

	def updateBatch(self, batch):
		"""
		updateBatch(batch)

		Updates the betweenness centralities after the batch `batch` of edge insertions.

		Parameters
		----------
		batch : list(networkit.dynamics.GraphEvent)
			Update the Betweenness based on a given list of graph events.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		(<_DynBetweenness*>(self._this)).updateBatch(_batch)

	def scores(self):
		""" 
		scores()

		Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		list(float)
			The betweenness scores calculated by run().
		"""
		return (<_DynBetweenness*>(self._this)).scores()

	def score(self, v):
		""" 
		score(v)

		Get the betweenness score of node `v` calculated by run().

		Parameters
		----------
		v : int
			A node.

		Returns
		-------
		float
			The betweenness score of node `v`.
		"""
		return (<_DynBetweenness*>(self._this)).score(v)

	def ranking(self):
		""" 
		ranking()
		
		Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score
		calculated by run().

		Returns
		-------
		list(tuple(int, float))
			A list of pairs.
		"""
		return (<_DynBetweenness*>(self._this)).ranking()

cdef extern from "<networkit/centrality/DynApproxBetweenness.hpp>":

	cdef cppclass _DynApproxBetweenness "NetworKit::DynApproxBetweenness"(_Algorithm):
		_DynApproxBetweenness(_Graph, double, double, bool_t, double) except +
		void update(_GraphEvent) except +
		void updateBatch(vector[_GraphEvent]) except +
		vector[double] scores() except +
		vector[pair[node, double]] ranking() except +
		double score(node) except +
		count getNumberOfSamples() except +

cdef class DynApproxBetweenness(Algorithm):
	""" 
	DynApproxBetweenness(G, epsilon=0.01, delta=0.1, storePredecessors=True, universalConstant=1.0)
	
	The algorithm approximates the betweenness of all vertices so that the scores are
	within an additive error epsilon with probability at least :math:`(1 - delta)`.
	The values are normalized by default.

	The algorithm approximates the betweenness of all vertices so that the scores are
	within an additive error epsilon with probability at least :math:`(1 - delta)`.
	The values are normalized by default.

	Parameters
	----------
	G : networkit.Graph
		The graph
	epsilon : float, optional
		Maximum additive error. Default: 0.01
	delta : float, optional
		Probability that the values are within the error guarantee. Default: 0.1
	storePredecessors : bool, optional
		Store lists of predecessors. Default: True
	universalConstant: float, optional
		The universal constant to be used in computing the sample size.
		It is 1 by default. Some references suggest using 0.5, but there
		is no guarantee in this case. Default: 1.0
	"""
	cdef Graph _G

	def __cinit__(self, Graph G, epsilon=0.01, delta=0.1, storePredecessors = True, universalConstant=1.0):
		self._G = G
		self._this = new _DynApproxBetweenness(G._this, epsilon, delta, storePredecessors, universalConstant)

	def update(self, ev):
		""" 
		update(ev)
		
		Updates the betweenness centralities after the edge insertions.

		Parameters
		----------
		ev : networkit.dynamics.GraphEvent
			Update the Betweenness based on a given graph event.
		"""
		(<_DynApproxBetweenness*>(self._this)).update(_GraphEvent(ev.type, ev.u, ev.v, ev.w))

	def updateBatch(self, batch):
		""" 
		updateBatch(batch)
		
		Updates the betweenness centralities after the batch `batch` of edge insertions.

		Parameters
		----------
		batch : list(networkit.dynamics.GraphEvent)
			Update the Betweenness based on a given list of graph events.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		(<_DynApproxBetweenness*>(self._this)).updateBatch(_batch)

	def scores(self):
		""" 
		scores()
		
		Get a vector containing the betweenness score for each node in the graph.

		Returns
		-------
		list(float)
			The betweenness scores calculated by run().
		"""
		return (<_DynApproxBetweenness*>(self._this)).scores()

	def score(self, v):
		""" 
		score(v)
		
		Get the betweenness score of node `v` calculated by run().

		Parameters
		----------
		v : int
			A node.

		Returns
		-------
		float
			The betweenness score of node `v`.
		"""
		return (<_DynApproxBetweenness*>(self._this)).score(v)

	def ranking(self):
		""" 
		ranking()
		
		Get a vector of pairs sorted into descending order. Each pair contains a node and the corresponding score
		calculated by run().

		Returns
		-------
		list(tuple(int, float))
			A list of pairs.
		"""
		return (<_DynApproxBetweenness*>(self._this)).ranking()

	def getNumberOfSamples(self):
		"""
		getNumberOfSamples()

		Get number of path samples used in last calculation.

		Returns
		-------
		int
			Number of samples.
		"""
		return (<_DynApproxBetweenness*>(self._this)).getNumberOfSamples()

cdef extern from "<networkit/centrality/DynBetweennessOneNode.hpp>":

	cdef cppclass _DynBetweennessOneNode "NetworKit::DynBetweennessOneNode":
		_DynBetweennessOneNode(_Graph, node) except +
		void run() nogil except +
		void update(_GraphEvent) except +
		void updateBatch(vector[_GraphEvent]) except +
		double getDistance(node, node) except +
		double getSigma(node, node) except +
		double getSigmax(node, node) except +
		double getbcx() except +

cdef class DynBetweennessOneNode:
	""" 
	DynBetweennessOneNode(G, x)
	
	Dynamic exact algorithm for updating the betweenness of a specific node.
	The algorithm aupdates the betweenness of a node after an edge insertions
	(faster than updating it for all nodes), based on the algorithm
	proposed by Bergamini et al. "Improving the betweenness centrality of a node by adding links"

	Parameters
	----------
	G : networkit.Graph
		The graph
	x : int
		The node for which you want to update betweenness
	"""
	cdef _DynBetweennessOneNode* _this
	cdef Graph _G

	def __cinit__(self, Graph G, node):
		self._G = G
		self._this = new _DynBetweennessOneNode(G._this, node)

	# this is necessary so that the C++ object gets properly garbage collected
	def __dealloc__(self):
		del self._this

	def run(self):
		with nogil:
			self._this.run()
		return self

	def update(self, ev):
		""" 
		update(ev)
		
		Updates the betweenness centralities after the edge insertions.

		Parameters
		----------
		ev : networkit.dynamics.GraphEvent
			Update the Betweenness based on a given graph event.
		"""
		self._this.update(_GraphEvent(ev.type, ev.u, ev.v, ev.w))

	def updateBatch(self, batch):
		""" 
		updateBatch(batch)
		
		Updates the betweenness centralities after the batch `batch` of edge insertions.

		Parameters
		----------
		batch : list(networkit.dynamics.GraphEvent)
			Update the Betweenness based on a given list of graph events.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		self._this.updateBatch(_batch)

	def getDistance(self, u, v):
		""" 
		getDistance(u, v)

		Returns the distance between node u and node v.

		Returns
		-------
		float
			Distance between u and v.
		"""
		return self._this.getDistance(u, v)

	def getSigma(self, u, v):
		"""
		getSigma(u, v)
		
		Returns the number of shortest paths between node u and node v.

		Returns
		-------
		int
			Number of shortest paths between u and v.
		"""
		return self._this.getSigma(u, v)

	def getSigmax(self, u, v):
		""" 
		getSigmax(u, v)
		
		Returns the number of shortest paths between node u and node v that go through x.
		
		Returns
		-------
		int
			Number of shortest paths between u and v that go through x.
		
		"""
		return self._this.getSigmax(u, v)

	def getbcx(self):
		""" 
		getbcx()
		
		Returns the betweenness centrality score of node x.

		Returns
		-------
		float
			Betweenness centrality score of x.
		"""
		return self._this.getbcx()

cdef extern from "<networkit/centrality/Closeness.hpp>" namespace "NetworKit":

	cdef enum _ClosenessVariant  "NetworKit::ClosenessVariant":
		standard = 0,
		generalized = 1

class ClosenessVariant:
	Standard = standard
	Generalized = generalized

cdef extern from "<networkit/centrality/Closeness.hpp>":

	cdef cppclass _Closeness "NetworKit::Closeness" (_Centrality):
		_Closeness(_Graph, bool, _ClosenessVariant) except +
		_Closeness(_Graph, bool, bool) except +

cdef class Closeness(Centrality):
	"""
	Closeness(G, normalized, checkConnectdedness)
	Closeness(G, normalized, variant)

	Constructs the Closeness class for the given Graph `G`. If the Closeness scores should not be normalized,
	set `normalized` to False. The run() method takes O(nm) time, where n is the number
	of nodes and m is the number of edges of the graph.

	Parameter :code:`variant` can be one of the following:

	- networkit.centrality.ClosenessVariant.Standard
	- networkit.centrality.ClosenessVariant.Generalized

	Parameters
	----------
	G : networkit.Graph
		The graph.
	normalized : bool
		Set this parameter to False if scores should not be normalized into an interval of [0,1].
		Normalization only works for unweighted graphs.
	checkConnectdedness : bool
		Set this parameter to True to also check if the graph is connected before computing closeness.
		Set this parameter to False to not check if the graph is connected (note: the standard definition
		of closeness works for connected graphs, choose this if the input graph is known to be connected).
	ClosenessVariant : networkit.centrality.ClosenessVariant
		Set this parameter to networkit.centrality.ClosenessVariant.Standard to use the standard
		definition of closeness, that is defined for connected graphs only; in this case, checkConnectdedness
		is automatically set to True.
		Set this parameter to networkit.centrality.ClosenessVariant.Generalized to use the generalized
		definition of closeness, that is defined for also non-connected graphs; in this case, checkConnectdedness
		is automatically set to False.
	"""

	def __cinit__(self, Graph G, normalized, third):
		self._G = G
		if isinstance(third, int):
			self._this = new _Closeness(G._this, normalized, <_ClosenessVariant> third)
		elif isinstance(third, bool):
			self._this = new _Closeness(G._this, normalized, <bool_t> third)
		else:
			raise Exception("Error: the third parameter must be either a bool or a ClosenessVariant")

cdef extern from "<networkit/centrality/ApproxCloseness.hpp>" namespace "NetworKit::ApproxCloseness":

	cdef enum _ClosenessType "NetworKit::ApproxCloseness::CLOSENESS_TYPE":
		INBOUND,
		OUTBOUND,
		SUM

class ClosenessType(object):
	Inbound = INBOUND
	Outbound = OUTBOUND
	Sum = SUM

cdef extern from "<networkit/centrality/ApproxCloseness.hpp>":

	cdef cppclass _ApproxCloseness "NetworKit::ApproxCloseness" (_Centrality):
		_ClosenessType type
		_ApproxCloseness(_Graph, count, float, bool_t, _ClosenessType type) except +
		vector[double] getSquareErrorEstimates() except +

cdef class ApproxCloseness(Centrality):
	""" 
	ApproxCloseness(G, nSamples, epsilon=0.1, normalized=False, type=networkit.centrality.ClosenessType.Outbound)

	Approximation of closeness centrality according to algorithm described in
  	Cohen et al., Computing Classic Closeness Centrality, at Scale.

	The algorithm approximates the closeness of all nodes in both directed and undirected graphs using a hybrid estimator.
	First, it takes nSamples samples. For these sampled nodes, the closeness is computed exactly. The pivot of each of the
	remaining nodes is the closest sampled node to it. If a node lies very close to its pivot, a sampling approach is used.
	Otherwise, a pivoting approach is used. Notice that the input graph has to be connected.

	Parameter :code:`type` can be one of the following:

	- networkit.centrality.ClosenessType.Inbound
	- networkit.centrality.ClosenessType.Outbound
	- networkit.centrality.ClosenessType.Sum

	Parameters
	----------
	G : networkit.Graph
		Input graph (undirected).
	nSamples : int
		User defined number of samples.
	epsilon : float, optional
		Parameter used for the error guarantee; it is also used to control when to use sampling and when to use pivoting. Default: 0.1
	normalized : bool, optional
		Normalize centrality values in interval [0,1]. Default: False
	type : networkit.centrality.ClosenessType, optional
		Use in- or outbound centrality or the sum of both (see paper) for computing closeness on directed graph. 
		If G is undirected, this can be ignored. Default: networkit.centrality.ClosenessType.Outbound
	"""

	#cdef _ApproxCloseness _this
	INBOUND = 0
	OUTBOUND = 1
	SUM = 2

	def __cinit__(self, Graph G, nSamples, epsilon=0.1, normalized=False, _ClosenessType type=OUTBOUND):
		self._G = G
		self._this = new _ApproxCloseness(G._this, nSamples, epsilon, normalized, type)

	def getSquareErrorEstimates(self):
		""" 
		getSquareErrorEstimates()
		
		Return a list containing the square error estimates for all nodes.

		Returns
		-------
		list(float)
			A list of square error estimate values.
		"""
		return (<_ApproxCloseness*>(self._this)).getSquareErrorEstimates()

cdef extern from "<networkit/centrality/DegreeCentrality.hpp>":

	cdef cppclass _DegreeCentrality "NetworKit::DegreeCentrality" (_Centrality):
		_DegreeCentrality(_Graph, bool_t normalized, bool_t outdeg, bool_t ignoreSelfLoops) except +

cdef class DegreeCentrality(Centrality):
	""" 
	DegreeCentrality(G, normalized=False, outDeg=True, ignoreSelfLoops=True)
	
	Node centrality index which ranks nodes by their degree.
	Optional normalization by maximum degree.  run() runs in O(n) time if ignoreSelfLoops is false or the graph 
	has no self-loops; otherwise it runs in O(m).

	Constructs the DegreeCentrality class for the given Graph `G`. If the scores should be normalized,
	then set `normalized` to True.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	normalized : bool, optional
		Normalize centrality values in the interval [0,1]. Default: False
	outdeg : bool, optional
		If set to true, computes the centrality based on out-degrees, otherwise based on the in-degrees. Default: True
	ignoreSelfLoops : bool, optional
		If set to true, self loops will not be taken into account. Default: True
	"""

	def __cinit__(self, Graph G, bool_t normalized=False, bool_t outDeg = True, bool_t ignoreSelfLoops=True):
		self._G = G
		self._this = new _DegreeCentrality(G._this, normalized, outDeg, ignoreSelfLoops)

cdef extern from "<networkit/centrality/HarmonicCloseness.hpp>":

	cdef cppclass _HarmonicCloseness "NetworKit::HarmonicCloseness" (_Centrality):
		_HarmonicCloseness(_Graph, bool_t) except +

cdef class HarmonicCloseness(Centrality):
	"""
	HarmonicCloseness(G, normalized=True)

	Constructs the HarmonicCloseness class for the given Graph `G`.
	If the harmonic closeness scores should not be normalized, set
	`normalized` to False.
	The run() method takes O(nm) time, where n is the number
	of nodes and m is the number of edges of the graph.

 	Parameters
 	----------
 	G : networkit.Graph
 		The input graph.
 	normalized : bool, optional
 		Set this parameter to False if scores should not be normalized 
		into an interval of [0,1]. Normalization only for unweighted graphs. Default: True
	"""

	def __cinit__(self, Graph G, normalized=True):
		self._G = G
		self._this = new _HarmonicCloseness(G._this, normalized)

cdef extern from "<networkit/centrality/TopCloseness.hpp>":

	cdef cppclass _TopCloseness "NetworKit::TopCloseness"(_Algorithm):
		_TopCloseness(_Graph G, count, bool_t, bool_t) except +
		node maximum() except +
		edgeweight maxSum() except +
		count iterations() except +
		count operations() except +
		vector[node] topkNodesList(bool_t) except +
		vector[edgeweight] topkScoresList(bool_t) except +


cdef class TopCloseness(Algorithm):
	"""
	TopCloseness(G, k=1, first_heu=True, sec_heu=True)
	
	Finds the top k nodes with highest closeness centrality faster than computing it for all nodes, 
	based on "Computing Top-k Closeness Centrality Faster in Unweighted Graphs", Bergamini et al., ALENEX16.
	The algorithms is based on two independent heuristics, described in the referenced paper. We recommend to 
	use :math:`first_heu = true` and :math:`second_heu = false` for complex networks and :math:`first_heu = true` 
	and :math:`second_heu = true` for street networks or networks with large diameters.

	Notes
	-----
	The worst case running time of the algorithm is :math:`O(nm)`, where n is the number of nodes and m is the number of edges.
	However, for most networks the empirical running time is :math:`O(m)`.

	Parameters
	----------
	G : networkit.Graph
		An unweighted graph.
	k : int
		Number of nodes with highest closeness that have to be found. 
		For example, if k = 10, the top 10 nodes with highest closeness will be computed.
	first_heu : bool, optional
		If true, the neighborhood-based lower bound is computed and nodes are sorted according to it. 
		If false, nodes are simply sorted by degree. Default: True
	sec_heu : bool, optional
		If true, the BFSbound is re-computed at each iteration. If false, BFScut is used. Default: True
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G, k=1, first_heu=True, sec_heu=True):
		self._G = G
		self._this = new _TopCloseness(G._this, k, first_heu, sec_heu)

	def topkNodesList(self, includeTrail=False):
		""" 
		topkNodesList(includeTrail=False)
		
		Returns: a list with the k nodes with highest closeness.
		WARNING: closeness centrality of some nodes below the top-k could be equal
	  	to the k-th closeness, we call them trail. Set the parameter includeTrail
	  	to true to also include those nodes but consider that the resulting vector
	  	could be longer than k.

		Parameters
		----------
		includeTrail : bool, optional
			Whether or not to include trail nodes. Default: False

		Returns
		-------
		list(int)
			The k nodes with highest closeness.
		"""
		return (<_TopCloseness*>(self._this)).topkNodesList(includeTrail)

	def topkScoresList(self, includeTrail=False):
		""" 
		topkScoresList(includeTrail=False)
		
		Returns: a list with the scores of the k nodes with highest closeness.
		WARNING: closeness centrality of some nodes below the top-k could be equal
  		to the k-th closeness, we call them trail. Set the parameter includeTrail
	  	to true to also include those centrality values but consider that the
	  	resulting vector could be longer than k.

		Parameters
		----------
		includeTrail : bool, optional
			Whether or not to include trail centrality value. Default: False

		Returns
		-------
		list(int)
			The k highest closeness scores.
		"""
		return (<_TopCloseness*>(self._this)).topkScoresList(includeTrail)

cdef extern from "<networkit/centrality/TopHarmonicCloseness.hpp>":

	cdef cppclass _TopHarmonicCloseness "NetworKit::TopHarmonicCloseness"(_Algorithm):
		_TopHarmonicCloseness(_Graph G, count, bool_t) except +
		vector[node] topkNodesList(bool_t) except +
		vector[edgeweight] topkScoresList(bool_t) except +


cdef class TopHarmonicCloseness(Algorithm):
	""" 
	TopHarmonicCloseness(G, k=1, useNBbound=True)

	Finds the top k nodes with highest harmonic closeness centrality faster
	than computing it for all nodes. The implementation is based on "Computing
	Top-k Centrality Faster in Unweighted Graphs", Bergamini et al., ALENEX16.
	The algorithm also works with weighted graphs but only if with the NBcut
	variation. We recommend to use useNBbound = False for complex (weighted)
	networks (or networks with small diameter) and useNBbound = True for
	unweighted street networks (or networks with large diameters). 
	
	Notes
	-----
	Notice that the worst case running time of the algorithm is O(nm), where n is the
	number of nodes and m is the number of edges. However, for most real-world
	networks the empirical running time is O(m).

	Parameters
	----------
	G : networkit.Graph
		The graph. If useNBbound is set to 'True', edge weights will be ignored.
	k : int, optional
		Number of nodes with highest closeness that have to be found. For example, if k = 10, the
		top 10 nodes with highest closeness will be computed. Default: 1
	useNBbound : bool, optional
		If True, the NBbound is re-computed at each iteration. If False, NBcut is used. The worst case 
		running time of the algorithm is :math:`O(nm)`, where n is the number of nodes and m is the number of edges.
		However, for most networks the empirical running time is :math:`O(m)`. Default: True
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G, k=1, useNBbound=False):
		self._G = G
		self._this = new _TopHarmonicCloseness(G._this, k, useNBbound)

	def topkNodesList(self, includeTrail=False):
		"""
		topkNodesList(includeTrail=False)

		Returns a list with the k nodes with highest harmonic closeness.
		WARNING: closeness centrality of some nodes below the top-k could be equal
		to the k-th closeness, we call them trail. Set the parameter includeTrail
		to true to also include those nodes but consider that the resulting vector
		could be longer than k.

		Parameters
		----------
		includeTrail : bool, optional
			Whether or not to include trail nodes. Default: False

		Returns
		-------
		list(int)
			The k nodes with highest harmonic closeness.
		"""
		return (<_TopHarmonicCloseness*>(self._this)).topkNodesList(includeTrail)

	def topkScoresList(self, includeTrail=False):
		"""
		topkScoresList(includeTrail=False)

		Returns a list with the scores of the k nodes with highest harmonic
		closeness. WARNING: closeness centrality of some nodes below the top-k
		could be equal to the k-th closeness, we call them trail. Set the
		parameter includeTrail to true to also include those centrality values
		but consider that the resulting vector could be longer than k.

		Parameters
		----------
		includeTrail : bool, optional
			Whether or not to include trail centrality value. Default: False

		Returns
		-------
		list(int)
			The k highest closeness harmonic scores.
		"""
		return (<_TopHarmonicCloseness*>(self._this)).topkScoresList(includeTrail)

cdef extern from "<networkit/centrality/DynTopHarmonicCloseness.hpp>":

	cdef cppclass _DynTopHarmonicCloseness "NetworKit::DynTopHarmonicCloseness"(_Algorithm):
		_DynTopHarmonicCloseness(_Graph G, count, bool_t) except +
		vector[pair[node, edgeweight]] ranking(bool_t) except +
		vector[node] topkNodesList(bool_t) except +
		vector[edgeweight] topkScoresList(bool_t) except +
		void update(_GraphEvent) except +
		void updateBatch(vector[_GraphEvent]) except +

cdef class DynTopHarmonicCloseness(Algorithm):
	""" 
	DynTopHarmonicCloseness(G, k=1, useBFSbound=True)
	
	Finds the top k nodes with highest harmonic closeness centrality faster
	than computing it for all nodes and updates them after a single or multiple
	edge update. The implementation is based on "Computing Top-k Closeness
	Centrality in Fully-dynamic Graphs", Bisenius et al., ALENEX18.
	The implementation is based on the static algorithms by Borassi et al.
	(complex networks) and Bergamini et al. (large-diameter networks).

	Notes
	-----
	The worst case running time of the algorithm is O(nm), where n is the number of nodes and m is the number of edges.
	However, for most networks the empirical running time is O(m).

	Parameters
	----------
	G : networkit.Graph
		An unweighted graph.
	k : int, optional
		Number of nodes with highest closeness that have to be found. 
		For example, if k = 10, the top 10 nodes with highest closeness will be computed. Default: 1
	useBFSbound : bool, optional
		If true, the BFSbound is re-computed at each iteration. If false, BFScut is used. Default: True
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G, k=1, useBFSbound=False):
		self._G = G
		self._this = new _DynTopHarmonicCloseness(G._this, k, useBFSbound)

	def ranking(self, includeTrail = False):
		""" 
		ranking(includeTrail = False)
		
		Returns: the ranking of the k most central nodes in the graph.
		WARNING: closeness centrality of some nodes below the top-k could be equal
		to the k-th closeness, we call them trail. Set the parameter includeTrail
		to true to also include those nodes but consider that the resulting vector
		could be longer than k.

		Parameters
		----------
		includeTrail : bool, optional
			Whether or not to include trail nodes. Default: False

		Returns
		-------
		int
			The ranking.
		"""
		return (<_DynTopHarmonicCloseness*>(self._this)).ranking(includeTrail)

	def topkNodesList(self, includeTrail = False):
		""" 
		topkNodesList(includeTrail = False)
		
		Returns: a list with the k nodes with highest harmonic closeness.
		WARNING: closeness centrality of some nodes below the top-k could be equal
		to the k-th closeness, we call them trail. Set the parameter includeTrail
		to true to also include those nodes but consider that the resulting vector
		could be longer than k.

		Parameters
		----------
		includeTrail : bool, optional
			Whether or not to include trail nodes. Default: False

		Returns
		-------
		list(int)
			The k nodes with highest harmonic closeness.
		"""
		return (<_DynTopHarmonicCloseness*>(self._this)).topkNodesList(includeTrail)

	def topkScoresList(self, includeTrail = False):
		""" 
		topkScoresList(includeTrail = False)
		
		Returns: a list with the scores of the k nodes with highest harmonic closeness.
		WARNING: closeness centrality of some nodes below the top-k could
		be equal to the k-th closeness, we call them trail. Set the parameter
		includeTrail to true to also include those centrality values but consider
		that the resulting vector could be longer than k.

		Parameters
		----------
		includeTrail : bool, optional
			Whether or not to include trail centrality value. Default: False

		Returns
		-------
		list(float)
			The k highest closeness harmonic scores.
		"""
		return (<_DynTopHarmonicCloseness*>(self._this)).topkScoresList(includeTrail)



	def update(self, ev):
		""" 
		update(ev)
		
		Updates the betweenness centralities after the edge insertions.

		Parameters
		----------
		ev : networkit.dynamics.GraphEvent
			Update the Betweenness based on a given graph event.
		"""
		(<_DynTopHarmonicCloseness*>(self._this)).update(_GraphEvent(ev.type, ev.u, ev.v, ev.w))

	def updateBatch(self, batch):
		""" 
		updateBatch(batch)
		
		Updates the betweenness centralities after the batch `batch` of edge insertions.

		Parameters
		----------
		batch : list(networkit.dynamics.GraphEvent)
			Update the Betweenness based on a given list of graph events.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		(<_DynTopHarmonicCloseness*>(self._this)).updateBatch(_batch)

cdef extern from "<networkit/centrality/LocalPartitionCoverage.hpp>":

	cdef cppclass _LocalPartitionCoverage "NetworKit::LocalPartitionCoverage" (_Centrality):
		_LocalPartitionCoverage(_Graph, _Partition) except +

cdef class LocalPartitionCoverage(Centrality):
	"""
	LocalPartitionCoverage(G, P)

	The local partition coverage is the amount of neighbors of a node u that are in the same partition as u.
	The running time of the run() method is O(m), where m is the number of edges in the graph.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	P : networkit.Partition
		The partition to use.
	"""
	cdef Partition _P

	def __cinit__(self, Graph G not None, Partition P not None):
		self._G = G
		self._P = P
		self._this = new _LocalPartitionCoverage(G._this, P._this)

cdef extern from "<networkit/centrality/GroupDegree.hpp>":

	cdef cppclass _GroupDegree "NetworKit::GroupDegree"(_Algorithm):
		_GroupDegree(_Graph G, count, bool_t) except +
		vector[node] groupMaxDegree() except +
		count getScore() except +
		count scoreOfGroup(vector[node]) except +

cdef class GroupDegree(Algorithm):
	"""
	GroupDegree(G, k = 1, countGroupNodes = True)

	Finds the group with the highest group degree centrality according to the
	definition proposed in 'The centrality of groups and classes' by Everett et
	al. (The Journal of mathematical sociology, 1999). This is a submodular but
	non monotone function so the algorithm can find a solution that is at least
	1/2 of the optimum. Worst-case running time is quadratic, but usually
	faster in real-world networks.
	The `countGroupNodes` option also count the nodes inside the group in the
	score, this make the group degree monotone and submodular and the algorithm
	is guaranteed to return a :math:`(1 - 1/e)`-approximation of the optimal solution.

	Parameters
	----------
		G : networkit.Graph
			The input graph.
		k: int, optional
			Size of the group of nodes. Default: 1
		countGroupNodes: bool, optional
			If nodes inside the group should be counted in the centrality score. Default: True
	"""
	cdef Graph _G

	def __cinit__(self, Graph G, k = 1, countGroupNodes = True):
		self._G = G
		self._this = new _GroupDegree(G._this, k, countGroupNodes)

	def groupMaxDegree(self):
		"""
		groupMaxDegree()

		Returns the group with maximum degree centrality.

		Returns
		-------
		list(int)
			The group of k nodes with highest degree centrality.
		"""
		return (<_GroupDegree*>(self._this)).groupMaxDegree()

	def getScore(self):
		"""
		getScore()

		Returns the score of the group with maximum degree centrality (i.e. the
		number of nodes outside the group that can be reached in one hop from at
		least one node in the group).

		Returns
		-------
		int
			The number of nodes outside the group that can be reached in one hop
			from at least one node in the group.
		"""
		return (<_GroupDegree*>(self._this)).getScore()

	def scoreOfGroup(self, vector[node] group):
		"""
		scoreOfGroup(vector[node] group)
		
		Returns the score of the given group.

		Parameters
		-----------
		group : list(int)
			List of nodes.

		Returns
		-------
		int
			The score of the given group.
		"""
		return (<_GroupDegree*>(self._this)).scoreOfGroup(group)

cdef extern from "<networkit/centrality/GedWalk.hpp>" namespace "NetworKit::GedWalk":

	cdef enum _BoundStrategy "NetworKit::GedWalk::BoundStrategy":
		no
		spectral
		geometric
		adaptiveGeometric

	cdef enum _GreedyStrategy "NetworKit::GedWalk::GreedyStrategy":
		lazy
		stochastic


class BoundStrategy(object):
	No = no
	Spectral = spectral
	Geometric = geometric
	AdaptiveGeometric = adaptiveGeometric

class GreedyStrategy(object):
	Lazy = lazy
	Stochastic = stochastic


cdef extern from "<networkit/centrality/GedWalk.hpp>":

	cdef cppclass _GedWalk "NetworKit::GedWalk"(_Algorithm):
		_GedWalk(_Graph G, count, double, double, _BoundStrategy, _GreedyStrategy, double) except +
		vector[node] groupMaxGedWalk() except +
		double getApproximateScore() except +
		double scoreOfGroup[InputIt](InputIt first, InputIt last, double epsilon) except +

cdef class GedWalk(Algorithm):
	"""
	GedWalk(Graph G, k = 1, epsilon = 0.1, alpha = -1.0, bs = networkit.centrality.BoundStrategy.Geometric, gs = GreedyStrategy.Lazy, spectralDelta = 0.5)

	Finds a group of `k` vertices with at least ((1 - 1/e) * opt - epsilon) GedWalk centrality
	score, where opt is the highest possible score. The algorithm is based on the paper "Group
	Centrality Maximization for Large-scale Graphs", Angriman et al., ALENEX20. It implements two
	independent greedy strategies (lazy and stochastic). Furthermore, it allows to compute the
	GedWalk score of a given set of nodes.

	Parameter :code:`bs` can be one of the following:

	- networkit.centrality.BoundStrategy.No
	- networkit.centrality.BoundStrategy.Spectral
	- networkit.centrality.BoundStrategy.Geometric
	- networkit.centrality.BoundStrategy.AdaptiveGeometric

	Parameter :code:`gs` can be one of the following:

	- networkit.centrality.GreedyStrategy.Lazy
	- networkit.centrality.GreedyStrategy.Stochastic

	Parameters
	----------
	G : networkit.Graph
		A (weakly) connected graph.
	k : int, optional
		The desired group size. Default: 1
	epsilon : float, optional
		Precision of the algorithm. Default: 0.1
	alpha : float, optional
		Exponent to compute the GedWalk score. Default: -1.0
	bs : networkit.centrality.BoundStrategy, optional
		Bound strategy to compute the GedWalk bounds.
		Default: networkit.centrality.BoundStrategy.Geometric
	gs : networkit.centrality.GreedyStrategy, optional
		Greedy strategy to be used (lazy or stochastic). 
		Default: networkit.centrality.GreedyStrategy.Lazy
	spectralDelta : float, optional
		Delta to be used for the spectral bound. Default: 0.5
	"""
	
	cdef Graph _G

	def __cinit__(self, Graph G, k = 1, epsilon = 0.1, alpha = -1.0, bs = BoundStrategy.Geometric,
			gs = GreedyStrategy.Lazy, spectralDelta = 0.5):

		self._G = G
		self._this = new _GedWalk(G._this, k, epsilon, alpha, bs, gs, spectralDelta)

	def __dealloc__(self):
		if self._this is not NULL:
			del self._this
			self._this = NULL

	def groupMaxGedWalk(self):
		"""
		groupMaxGedWalk()

		Returns the computed group.

		Returns
		-------
		list(int)
			The computed group.
		"""
		return (<_GedWalk*>(self._this)).groupMaxGedWalk()

	def getApproximateScore(self):
		"""
		getApproximateScore()

		Returns the GedWalk score of the computed group.

		Returns
		-------
		float
			The GedWalk score of the computed group.
		"""
		return (<_GedWalk*>(self._this)).getApproximateScore()

	def scoreOfGroup(self, group, epsilon = 0.1):
		"""
		scoreOfGroup(group, epsilon = 0.1)
		
		Returns the GedWalk score of the input group.

		Parameters
		----------
		group : list(int)
			The input group.
		epsilon : float, optional
			The precision of the score to be computed. Default: 0.1

		Returns
		-------
		float
			An epsilon-approximation of the GedWalk score of the input group.
		"""
		cdef vector[node] groupVec

		try:
			groupVec = <vector[node]?>group
		except TypeError:
			raise RuntimeError("Error, group must be a list of nodes.")
		return (<_GedWalk*>(self._this)).scoreOfGroup[vector[node].iterator](groupVec.begin(), groupVec.end(), epsilon)


cdef extern from "<networkit/centrality/ApproxGroupBetweenness.hpp>":

	cdef cppclass _ApproxGroupBetweenness "NetworKit::ApproxGroupBetweenness" (_Algorithm):
		_ApproxGroupBetweenness(_Graph, count, double) except +
		vector[node] groupMaxBetweenness() except +
		count scoreOfGroup(vector[node]) except +

cdef class ApproxGroupBetweenness(Algorithm):
	"""
	ApproxGroupBetweenness(G, groupSize, epsilon)

	Constructs the ApproxGroupBetweenness class for a given undirected Graph
	`G`.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	groupSize : int
		The desired size of the group.
	epsilon : float
		Determines the accuracy of the approximation.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G, groupSize, epsilon):
		self._G = G
		self._this = new _ApproxGroupBetweenness(G._this, groupSize, epsilon)

	def groupMaxBetweenness(self):
		"""
		groupMaxBetweenness()
		
		Get a vector of nodes containing the set of nodes with apporoximated
		maximum group betweenness.

		Returns
		-------
		list(int)
			The group of nodes with highest approximated group betweenness.
		"""
		return (<_ApproxGroupBetweenness*>(self._this)).groupMaxBetweenness()

	def scoreOfGroup(self, vector[node] group):
		"""
		scoreOfGroup(group)
		
		Returns the score of the given group.

		Parameters
		----------
		group : list(int)
			Set of nodes.

		Returns
		-------
		int
			The score of the given group.
		"""
		return (<_ApproxGroupBetweenness*>(self._this)).scoreOfGroup(group)

cdef extern from "<networkit/centrality/GroupCloseness.hpp>":

		cdef cppclass _GroupCloseness "NetworKit::GroupCloseness"(_Algorithm):
			_GroupCloseness(_Graph G, count, count) except +
			vector[node] groupMaxCloseness() except +
			double computeFarness(vector[node], count) except +
			double scoreOfGroup(vector[node]) except +

cdef class GroupCloseness(Algorithm):
	"""
	GroupCloseness(G, k=1, H=0)
	
	Finds the group of nodes with highest (group) closeness centrality. The algorithm is the one 
	proposed in Bergamini et al., ALENEX 2018 and finds a solution that is a (1-1/e)-approximation of the optimum.
	The worst-case running time of this approach is quadratic, but usually much faster in practice.

	Parameters
	----------
	G : networkit.Graph
		An unweighted graph.
	k : int, optional
		Size of the group. Default: 1
	H : int, optional
		If equal 0, simply runs the algorithm proposed in Bergamini et al.. 
		If > 0, interrupts all BFSs after H iterations (suggested for very large networks). Default: 0
	"""
	cdef Graph _G

	def __cinit__(self,  Graph G, k=1, H=0):
		self._G = G
		self._this = new _GroupCloseness(G._this, k, H)

	def groupMaxCloseness(self):
		"""
		groupMaxCloseness()
		
		Returns the group with maximum closeness centrality.

		Returns
		-------
		list(int)
			The group of k nodes with maximum closeness centrality.
		"""
		return (<_GroupCloseness*>(self._this)).groupMaxCloseness()

	def computeFarness(self, S, H=0):
		""" 
		computeFarness(S, H=0)
		
		Computes farness (i.e., inverse of the closeness) for a given group (stopping after H iterations if H > 0).

		Parameters
		----------
		S : list(int)
			Group to compute farness on.
		H : int, optional
			If equal 0, simply runs the algorithm proposed in Bergamini et al.. 
			If > 0, interrupts after H iterations (suggested for very large networks). Default: 0

		Returns
		-------
		float
			Farness value for node group.
		"""
		return (<_GroupCloseness*>(self._this)).computeFarness(S, H)

	def scoreOfGroup(self, group):
		"""
		scoreOfGroup(group)
		
		Computes the group closeness score of the given group.

		Parameters
		----------
		group : list(int)
			List of nodes.

		Returns
		-------
		float
			The group closeness score of the given group.
		"""
		return (<_GroupCloseness*>(self._this)).scoreOfGroup(group)


cdef extern from "<networkit/centrality/GroupClosenessGrowShrink.hpp>":
	cdef cppclass _GroupClosenessGrowShrink "NetworKit::GroupClosenessGrowShrink"(_Algorithm):
		_GroupClosenessGrowShrink(_Graph, const vector[node] &group, bool_t extended, count insertions) except +
		vector[node] groupMaxCloseness() except +
		count numberOfIterations() except +
		count maxIterations

cdef class GroupClosenessGrowShrink(Algorithm):
	"""
	GroupClosenessGrowShrink(Graph G, group, extended = False, insertions = 0)

	Finds a group of nodes with high group closeness centrality. This is the Grow-Shrink
	algorithm presented in Angriman et al. "Local Search for Group Closeness Maximization on Big
	Graphs" IEEE BigData 2019. The algorithm takes as input a connected, unweighted, undirected
	graph and an arbitrary group of nodes, and improves the group closeness of the given
	group by performing vertex exchanges.

	Parameters
	----------
	G : networkit.Graph
		A connected, undirected, unweighted graph.
	group : list(int)
		The initial group of nodes.
	extended : bool, optional
		Set this parameter to true for the Extended Grow-Shrink algorithm (i.e.,
		swaps are not restricted to only neighbors of the group). Default: False
	insertions : int, optional
		Number of consecutive node insertions and removal per iteration. Let this
		parameter to zero to use Diameter(G)/sqrt(k) nodes (where k is the size of the group). Default: 0
	"""
	cdef Graph _G

	def __cinit__(self, Graph G, group, extended = False, insertions = 0):
		self._G = G
		self._this = new _GroupClosenessGrowShrink(G._this, group, extended, insertions)

	def groupMaxCloseness(self):
		"""
		groupMaxCloseness()

		Returns the computed group.

		Returns
		-------
		list(int)
			The computed group.
		"""
		return (<_GroupClosenessGrowShrink*>(self._this)).groupMaxCloseness()

	def numberOfIterations(self):
		"""
		Return the total number of iterations performed by the algorithm.

		Returns
		-------
		int
			Total number of iterations performed by the algorithm.
		"""
		return (<_GroupClosenessGrowShrink*>(self._this)).numberOfIterations()


cdef extern from "<networkit/centrality/GroupClosenessLocalSwaps.hpp>":
	cdef cppclass _GroupClosenessLocalSwaps "NetworKit::GroupClosenessLocalSwaps"(_Algorithm):
		_GroupClosenessLocalSwaps(_Graph, const vector[node] &group, count maxSwaps) except +
		vector[node] groupMaxCloseness() except +
		count numberOfSwaps() except +

cdef class GroupClosenessLocalSwaps(Algorithm):
	"""
	GroupClosenessLocalSwaps(Graph G, group, maxSwaps = 0)

	Finds a group of nodes with high group closeness centrality. This is
	the LS-restrict algorithm presented in Angriman et al. "Local Search
	for Group Closeness Maximization on Big Graphs" IEEE BigData 2019. The
	algorithm takes as input a graph and an arbitrary group of nodes, and
	improves the group closeness of the given
	group by performing vertex exchanges.

	Parameters
	----------
	G : networkit.Graph
		An undirected, unweighted graph.
	group : list(int)
		The initial group of nodes.
	maxSwaps : int, optional
		Maximum number of vertex exchanges allowed. Default: 0
	"""
	cdef Graph _G

	def __cinit__(self, Graph G, group, maxSwaps = 0):
		self._G = G
		self._this = new _GroupClosenessLocalSwaps(G._this, group, maxSwaps)

	def groupMaxCloseness(self):
		"""
		groupMaxCloseness()

		Returns the computed group.

		Returns
		-------
		list(int)
			The computed group.
		"""
		return (<_GroupClosenessLocalSwaps*>(self._this)).groupMaxCloseness()

	def numberOfSwaps(self):
		"""
		numberOfSwaps()

		Return the total number of vertex exchanges performed by the algorithm.

		Returns
		-------
		int
			Total number of vertex exchanges performed by the algorithm.
		"""
		return (<_GroupClosenessLocalSwaps*>(self._this)).numberOfSwaps()


cdef extern from "<networkit/centrality/GroupHarmonicCloseness.hpp>":
	cdef cppclass _GroupHarmonicCloseness "NetworKit::GroupHarmonicCloseness"(_Algorithm):
		_GroupHarmonicCloseness(_Graph G, count) except +
		vector[node] groupMaxHarmonicCloseness() except +
		@staticmethod
		double scoreOfGroup[InputIt](_Graph G, InputIt first, InputIt last) except +

cdef class GroupHarmonicCloseness(Algorithm):
	"""
	GroupHarmonicCloseness(Graph G, k = 1)

	Approximation algorithm for the group-harmonic maximization problem. The
	computed solutions have a guaranteed :math:`\lambda(1 - \\frac{1}{2e})`
	(directed graphs) and :math:`\lambda(1 - \\frac{1}{e})/2` (undirected graphs)
	approximation ratio, where :math:`\lambda` is the ratio between the minimal and
	the maximal edge weight. The algorithm is the one proposed in Angriman et
	al., ALENEX 2021. 
	
	Notes
	-----
	The worst-case running time of this approach is
	quadratic, but usually much faster in practice.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	k : int, optional
		Size of the group of nodes. Default: 1
	"""
	cdef Graph _G

	def __cinit__(self, Graph G, k = 1):
		self._G = G
		self._this = new _GroupHarmonicCloseness(G._this, k)

	def groupMaxHarmonicCloseness(self):
		"""
		groupMaxHarmonicCloseness()

		Returns the computed group.

		Returns
		-------
		list(int)
			The computed group.
		"""
		return (<_GroupHarmonicCloseness*>(self._this)).groupMaxHarmonicCloseness()

	@staticmethod
	def scoreOfGroup(Graph graph, vector[node] inputGroup):
		"""
		scoreOfGroup(Graph graph, vector[node] inputGroup)
		
		Computes the group-harmonic score of the input group.

		Parameters
		----------
		graph : networkit.Graph
			The input graph.
		inputGroup : list(int)
			The input group of nodes.

		Returns
		-------
		float
			The group-harmonic score of the input group.
		"""
		return _GroupHarmonicCloseness.scoreOfGroup[vector[node].iterator](
				graph._this, inputGroup.begin(), inputGroup.end())

cdef extern from "<networkit/centrality/GroupClosenessLocalSearch.hpp>":
	cdef cppclass _GroupClosenessLocalSearch "NetworKit::GroupClosenessLocalSearch"(_Algorithm):
		_GroupClosenessLocalSearch(_Graph, const vector[node] &group, bool_t runGrowShrink, count maxIterations) except +
		_GroupClosenessLocalSearch(_Graph, const vector[node] &group, bool_t runGrowShrink) except +
		_GroupClosenessLocalSearch(_Graph, const vector[node] &group) except +
		vector[node] groupMaxCloseness() except +
		count numberOfIterations() except +
		count maxIterations

cdef class GroupClosenessLocalSearch(Algorithm):
	"""
	GroupClosenessLocalSearch(Graph G, group, runGrowShrink, maxIterations)

	Local search approximation algorithm for Group Closeness Maximization presented in
	"Group-Harmonic and Group-Closeness Maximization – Approximation and Engineering", Angriman
	et al., ALENEX 2021. The algorithm evaluates all possible swaps between a vertex in the group
	and the vertices outside, and performs a swap iff the decrement in farness is at least :math:`(1 -
	1 / (k \\cdot (n - k)))`, where `k` is the number of vertices in the group. Thus,
	even in a best-case scenario the time complexity of this algorithm is :math:`O(n \\cdot k)`. To
	keep the number of swaps low, it is recommended to use this algorithm as a refinement step of
	an already good solution computed by a faster algorithm e.g., greedy (GroupCloseness), or
	GrowShrink (GroupClosenessGrowShrink). In undirected graphs the approximation ratio is 1/5,
	on directed graphs it has not been demonstrated.

	Parameters
	----------
	G : networkit.Graph
		A connected, undirected, unweighted graph.
	group : list(int)
		The initial group of nodes.
	useGrowShrink : bool
		Whether or not to run the GrowShrink algorithm on the initial group.
	maxIterations : int
		Maximum number of swaps allowed. Prevents the algorithm from performing
		too many swaps by giving up the approximation guarantee.
	"""
	cdef Graph _G

	def __cinit__(self, Graph G not None, group, runGrowShrink, maxIterations):
		self._G = G
		self._this = new _GroupClosenessLocalSearch(G._this, group, runGrowShrink, maxIterations)

	def __cinit__(self, Graph G not None, group, runGrowShrink):
		self._G = G
		self._this = new _GroupClosenessLocalSearch(G._this, group, runGrowShrink)

	def __cinit__(self, Graph G not None, group):
		self._G = G
		self._this = new _GroupClosenessLocalSearch(G._this, group)

	def groupMaxCloseness(self):
		"""
		groupMaxCloseness()

		Returns the computed group.

		Returns
		-------
		list(int)
			The computed group.
		"""
		return (<_GroupClosenessLocalSearch*>(self._this)).groupMaxCloseness()

	def numberOfIterations(self):
		"""
		numberOfIterations()

		Return the total number of iterations performed by the algorithm.

		Returns
		-------
		int
			Total number of iterations performed by the algorithm.
		"""
		return (<_GroupClosenessLocalSearch*>(self._this)).numberOfIterations()


cdef extern from "<networkit/centrality/KPathCentrality.hpp>":

	cdef cppclass _KPathCentrality "NetworKit::KPathCentrality" (_Centrality):
		_KPathCentrality(_Graph, double, count) except +

cdef class KPathCentrality(Centrality):
	"""
	KPathCentrality(G, alpha=0.2, k=0)

	Constructs the K-Path Centrality class for the given Graph `G`.

 	Parameters
 	----------
 	G : networkit.Graph
 		The graph.
 	alpha : float, optional
		Tradeoff between runtime and precision with -0.5: maximum precision, maximum runtime
 		and 0.5: lowest precision, lowest runtime. Default: 0.2
	k: int, optional
		maximum length of paths. Default: 0
	"""

	def __cinit__(self, Graph G, alpha=0.2, k=0):
		self._G = G
		self._this = new _KPathCentrality(G._this, alpha, k)

cdef extern from "<networkit/centrality/KatzCentrality.hpp>" namespace "NetworKit":

	cdef enum _EdgeDirection"NetworKit::EdgeDirection":
		OutEdges = 0
		InEdges = 1

class EdgeDirection(object):
	inEdges = InEdges
	outEdges = OutEdges

cdef extern from "<networkit/centrality/KatzCentrality.hpp>":

	cdef cppclass _KatzCentrality "NetworKit::KatzCentrality" (_Centrality):
		_KatzCentrality(_Graph, double, double, double) except +
		_EdgeDirection edgeDirection

cdef class KatzCentrality(Centrality):
	"""
	KatzCentrality(G, alpha=0, beta=0.1, tol=1e-8)

	Constructs a KatzCentrality object for the given Graph `G`.
	Each iteration of the algorithm requires O(m) time.
	The number of iterations depends on how long it takes to reach the convergence
	(and therefore on the desired tolerance `tol`).

 	Parameters
 	----------
 	G : networkit.Graph
 		The input graph.
 	alpha : float, optional
		Damping of the matrix vector product result, must be non negative.
		Leave this parameter to 0 to use the default value :math:`1 / (max_degree + 1)`. Default: 0
	beta : float, optional
		Constant value added to the centrality of each vertex. Default: 0.1
	tol : float, optional
		The tolerance for convergence. Default: 1e-8
	"""

	def __cinit__(self, Graph G, alpha=0, beta=0.1, tol=1e-8):
		self._G = G
		self._this = new _KatzCentrality(G._this, alpha, beta, tol)

	property edgeDirection:
		"""
		Property :code:`edgeDirection` can be one of the following:

		- networkit.centrality.EdgeDirection.inEdges
		- networkit.centrality.EdgeDirection.outEdges

		Default: networkit.centrality.EdgeDirection.inEdges
		"""
		def __get__(self):
			""" Get the used edge direction. """
			return (<_KatzCentrality*>(self._this)).edgeDirection
		def __set__(self, _EdgeDirection edgeDirection):
			""" Use a different edge direction. """
			(<_KatzCentrality*>(self._this)).edgeDirection = edgeDirection

cdef extern from "<networkit/centrality/DynKatzCentrality.hpp>":

	cdef cppclass _DynKatzCentrality "NetworKit::DynKatzCentrality" (_Centrality):
		_DynKatzCentrality(_Graph G, count, bool_t, double) except +
		void update(_GraphEvent) except +
		void updateBatch(vector[_GraphEvent]) except +
		node top(count) except +
		double bound(node) except +
		bool_t areDistinguished(node, node) except +

cdef class DynKatzCentrality(Centrality):
	""""
	DynKatzCentrality(G, k, groupOnly=False, tolerance=1e-9)

	Finds the top-k nodes with highest Katz centrality.

 	Parameters
 	----------
 	G : networkit.Graph
 		The input graph.
 	k : float
		The number k for which we want to find the top-k nodes with highest Katz centrality.
	groupOnly : bool, optional
		Set whether the update will only update top-k nodes. Default: False
	tolerance : float, optional
		The tolerance for convergence. Default: 1e-9
	"""

	def __cinit__(self, Graph G, k, groupOnly=False, tolerance=1e-9):
		self._G = G
		self._this = new _DynKatzCentrality(G._this, k, groupOnly, tolerance)

	def update(self, ev):
		""" 
		update(ev)
		
		Updates the centralities after the edge insertions.

		Parameters
		----------
		ev : networkit.dynamics.GraphEvent
			Update the Betweenness based on a given graph event.
		"""
		(<_DynKatzCentrality*>(self._this)).update(_GraphEvent(ev.type, ev.u, ev.v, ev.w))

	def updateBatch(self, batch):
		""" 
		updateBatch(batch)
		
		Updates the centralities after the batch `batch` of edge insertions.

		Parameters
		----------
		batch : list(networkit.dynamics.GraphEvent)
			Update the Betweenness based on a given list of graph events.
		"""
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		(<_DynKatzCentrality*>(self._this)).updateBatch(_batch)

	def top(self, n=0):
		""" 
		top(n=0)
		
		Returns the top n nodes.

		Parameters
		----------
		n : int, optional
			If set, retrieve n top-nodes. If not set, all nodes are retrieved. Default: 0

		Returns
		-------
		list(int)
			List of nodes with top-n centrality-scores.
		"""
		return (<_DynKatzCentrality*>(self._this)).top(n)

	def bound(self, v):
		""" 
		bound(v)
		
		Returns the (upper) bound of the centrality of node `v`.

		Parameters
		----------
		v : int
			Node in the graph.

		Returns
		-------
		float
			Upper bound of node `v`.
		"""		
		return (<_DynKatzCentrality*>(self._this)).bound(v)

	def areDistinguished(self, u, v):
		""" 
		areDistinguished(u, v)
		
		Returns true if the bounds are sharp enough to rank two nodes against each other.

		Parameters
		----------
		u : int
			Node in the graph.
		v : int
			Node in the graph.

		Returns
		-------
		bool
			True if the bounds are sharp enough to rank two nodes against each other.
		"""	
		return (<_DynKatzCentrality*>(self._this)).areDistinguished(u, v)

cdef extern from "<networkit/centrality/LocalClusteringCoefficient.hpp>":

	cdef cppclass _LocalClusteringCoefficient "NetworKit::LocalClusteringCoefficient" (_Centrality):
		_LocalClusteringCoefficient(_Graph, bool_t) except +

cdef class LocalClusteringCoefficient(Centrality):
	"""
	LocalClusteringCoefficient(G, turbo=False)

	Constructs the LocalClusteringCoefficient class for the given Graph `G`. If the local clustering coefficient values should be normalized,
	then set `normalized` to True. The graph may not contain self-loops.

	There are two algorithms available. The trivial (parallel) algorithm needs only a small amount of additional memory.
	The turbo mode adds a (sequential, but fast) pre-processing step using ideas from Triangle Listing Algorithms: Back 
	from the Diversion (Mark Ortmann and Ulrik Brandes). This reduces the running time significantly for most graphs. 
	However, the turbo mode needs O(m) additional memory. In practice this should be a bit less than half of the memory 
	that is needed for the graph itself. The turbo mode is particularly effective for graphs with nodes of very high degree 
	and a very skewed degree distribution.

 	Parameters
 	----------
 	G : networkit.Graph
 		The input graph.
	turbo : bool, optional
		If the turbo mode shall be activated. Default: False
	"""

	def __cinit__(self, Graph G, bool_t turbo = False):
		self._G = G
		self._this = new _LocalClusteringCoefficient(G._this, turbo)

cdef extern from "<networkit/centrality/LocalSquareClusteringCoefficient.hpp>":

	cdef cppclass _LocalSquareClusteringCoefficient "NetworKit::LocalSquareClusteringCoefficient" (_Centrality):
		_LocalSquareClusteringCoefficient(_Graph) except +


cdef class LocalSquareClusteringCoefficient(Centrality):
	"""
	LocalSquareClusteringCoefficient(G)

	Constructs the LocalSquareClusteringCoefficient class for the given Graph `G`.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _LocalSquareClusteringCoefficient(G._this)


cdef extern from "<networkit/centrality/Sfigality.hpp>":

	cdef cppclass _Sfigality "NetworKit::Sfigality" (_Centrality):
		_Sfigality(_Graph) except +

cdef class Sfigality(Centrality):
	"""
	Sfigality(G)

	Sfigality is a new type of node centrality measures that is high if neighboring nodes have a higher degree, e.g. in social networks, 
	if your friends have more friends than you. Formally this means: :math:`\sigma(u) = \\frac{| \{ v: \{u,v\} \in E, deg(u) < deg(v) \} |}{ deg(u) }`

 	Parameters
 	----------
 	G : networkit.Graph
 		The input graph.
	"""

	def __cinit__(self, Graph G):
		self._G = G
		self._this = new _Sfigality(G._this)

cdef extern from "<networkit/centrality/PermanenceCentrality.hpp>":

	cdef cppclass _PermanenceCentrality "NetworKit::PermanenceCentrality"(_Algorithm):
		_PermanenceCentrality(const _Graph& G, const _Partition& P) except +
		double getIntraClustering(node u) except +
		double getPermanence(node u) except +

cdef class PermanenceCentrality(Algorithm):
	"""
	PermanenceCentrality(Graph G, Partition P)

	This centrality measure measure how well a vertex belongs to its community. The values are calculated on the fly, the 
	partion may be changed in between the requests.
	For details see	On the permanence of vertices in network communities (Tanmoy Chakraborty, Sriram Srinivasan, 
	Niloy Ganguly, Animesh Mukherjee, and Sanjukta Bhowmick. 2014)

 	Parameters
 	----------
 	G : networkit.Graph
 		The input graph.
	P : networkit.Partition
		Partition for graph G.
	"""
	cdef Graph _G
	cdef Partition _P

	def __cinit__(self, Graph G, Partition P):
		self._this = new _PermanenceCentrality(G._this, P._this)
		self._G = G
		self._P = P

	def getIntraClustering(self, node u):
		""" 
		getIntraClustering(u)
		
		Returns intra clustering for node u.

		Parameters
		----------
		u : int
			Node in the graph.

		Returns
		-------
		float
			Intra clustering value for node u.
		"""
		return (<_PermanenceCentrality*>(self._this)).getIntraClustering(u)

	def getPermanence(self, node u):
		""" 
		getPermanence(u)
		
		Returns permanence centrality for node u.

		Parameters
		----------
		u : int
			Node in the graph.

		Returns
		-------
		float
			Permanence centrality value for node u.
		"""
		return (<_PermanenceCentrality*>(self._this)).getPermanence(u)


cdef extern from "<networkit/centrality/LaplacianCentrality.hpp>":

	cdef cppclass _LaplacianCentrality "NetworKit::LaplacianCentrality" (_Centrality):
		_LaplacianCentrality(_Graph, bool_t) except +

cdef class LaplacianCentrality(Centrality):
	""" 
	LaplacianCentrality(G, normalized=False)
	
	Computes the Laplacian centrality of the graph.

	Notes
	-----
	The implementation is a simplification of the original algorithm proposed by Qi et al. in
	"Laplacian centrality: A new centrality measure for weighted networks".

	See https://dl.acm.org/citation.cfm?id=2181343.2181780 for details.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	normalized : bool, optional
		Whether scores should be normalized by the energy of the full graph. Default: False
	"""

	def __cinit__(self, Graph G, normalized = False):
		self._G = G
		self._this = new _LaplacianCentrality(G._this, normalized)

cdef extern from "<networkit/centrality/CoreDecomposition.hpp>":

	cdef cppclass _CoreDecomposition "NetworKit::CoreDecomposition" (_Centrality):
		_CoreDecomposition(_Graph, bool_t, bool_t, bool_t) except +
		_Cover getCover() except +
		_Partition getPartition() except +
		index maxCoreNumber() except +
		vector[node] getNodeOrder() except +

cdef class CoreDecomposition(Centrality):
	""" 
	CoreDecomposition(G, normalized=False, enforceBucketQueueAlgorithm=False, storeNodeOrder = False)

	Computes k-core decomposition of a graph. The graph may not contain self-loops.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	normalized : bool, optional
		Divide each core number by the maximum degree. Default: False
	enforceBucketQueueAlgorithm : bool, optional
		Enforce switch to sequential algorithm. Default: False
	storeNodeOrder : bool, optional
		If set to True, the order of the nodes in ascending order of the cores is stored and can later be returned using getNodeOrder(). 
		Enforces the sequential bucket priority queue algorithm. Default: False

	"""

	def __cinit__(self, Graph G, bool_t normalized=False, bool_t enforceBucketQueueAlgorithm=False, bool_t storeNodeOrder = False):
		self._G = G
		self._this = new _CoreDecomposition(G._this, normalized, enforceBucketQueueAlgorithm, storeNodeOrder)

	def maxCoreNumber(self):
		""" 
		maxCoreNumber()

		Get maximum core number.

		Returns
		-------
		int
			The maximum core number.
		"""
		return (<_CoreDecomposition*>(self._this)).maxCoreNumber()

	def getCover(self):
		""" 
		getCover()
		
		Get the k-cores as cover.

		Returns
		-------
		list(int)
			The k-cores as sets of nodes, indexed by k.
		"""
		return Cover().setThis((<_CoreDecomposition*>(self._this)).getCover())

	def getPartition(self):
		""" 
		getPartition()
		
		Get the k-shells as a partition object.

		Returns
		-------
		networkit.Partition
			The k-shells.
		"""
		return Partition().setThis((<_CoreDecomposition*>(self._this)).getPartition())

	def getNodeOrder(self):
		"""
		getNodeOrder()

		Get the node order. This is only possible when storeNodeOrder was set.

		Returns
		-------
		list(int)
			The nodes sorted by increasing core number.
		"""
		return (<_CoreDecomposition*>(self._this)).getNodeOrder()

cdef extern from "<networkit/centrality/EigenvectorCentrality.hpp>":

	cdef cppclass _EigenvectorCentrality "NetworKit::EigenvectorCentrality" (_Centrality):
		_EigenvectorCentrality(_Graph, double tol) except +

cdef class EigenvectorCentrality(Centrality):
	"""	
	EigenvectorCentrality(G, tol=1e-9)
	
	Computes the leading eigenvector of the graph's adjacency matrix (normalized in 2-norm).
	Interpreted as eigenvector centrality score.

	Constructs the EigenvectorCentrality class for the given Graph `G`. `tol` defines the tolerance for convergence.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	tol : float, optional
		The tolerance for convergence.
	"""

	def __cinit__(self, Graph G, double tol=1e-9):
		self._G = G
		self._this = new _EigenvectorCentrality(G._this, tol)

cdef extern from "<networkit/centrality/PageRank.hpp>" namespace "NetworKit::PageRank":

	#cdef enum Norm:
	cdef enum _Norm "NetworKit::PageRank::Norm":
		L1Norm = 0
		L2Norm = 1

	cdef enum _SinkHandling "NetworKit::PageRank::SinkHandling":
		NO_SINK_HANDLING
		DISTRIBUTE_SINKS

class Norm(object):
	l1norm = L1Norm
	l2norm = L2Norm

class SinkHandling(object):
	NoSinkHandling = NO_SINK_HANDLING
	DistributeSinks = DISTRIBUTE_SINKS


cdef extern from "<networkit/centrality/PageRank.hpp>":

	cdef cppclass _PageRank "NetworKit::PageRank" (_Centrality):
		_PageRank(_Graph, double damp, double tol, bool_t normalized, _SinkHandling distributeSinks) except +
		count numberOfIterations() except +
		_Norm norm
		count maxIterations

cdef class PageRank(Centrality):
	""" 
	PageRank(G, damp=0.85, tol=1e-9, normalized=False, distributeSinks=SinkHandling.NoSinkHandling)
	
	Compute PageRank as node centrality measure. In the default mode this computation is in line
 	with the original paper "The PageRank citation ranking: Bringing order to the web." by L. Brin et al (1999).
 	In later publications ("PageRank revisited." by M. Brinkmeyer et al. (2005) amongst others) sink-node handling
 	was added for directed graphs in order to comply with the theoretical assumptions by the underlying 
 	Markov chain model. This can be activated by setting the matching parameter to true. By default 
 	this is disabled, since it is an addition to the original definition.
	
 	Page-Rank values can also be normalized by post-processed according to "Comparing Apples and
 	Oranges: Normalized PageRank for Evolving Graphs" by Berberich et al. (2007). This decouples
 	the PageRank values from the size of the input graph. To enable this, set the matching parameter
 	to true. Note that, sink-node handling is automatically activated if normalization is used.

	Parameter :code:`distributeSinks` can be one of the following:

	- networkit.centrality.SinkHandling.NoSinkHandling
	- networkit.centrality.SinkHandling.DistributeSinks

	Parameters
	----------
	G : networkit.Graph
		Graph to be processed.
	damp : float, optional
		Damping factor of the PageRank algorithm. Default: 0.85
	tol : float, optional
		Error tolerance for PageRank iteration. Default: 1e-8
	distributeSinks: networkit.centrality.SinkHandling, optional
		Set to distribute PageRank values for sink nodes. Default: SinkHandling.NoSinkHandling
	normalized : bool, optional
		If the results should be normalized by the lower bound of scores. 
		This decouples the PageRank values from the size of the input graph. Default: False
	"""

	def __cinit__(self, Graph G, double damp=0.85, double tol=1e-8, bool_t normalized=False, distributeSinks=SinkHandling.NoSinkHandling):
		self._G = G
		self._this = new _PageRank(G._this, damp, tol, normalized, distributeSinks)

	def numberOfIterations(self):
		"""
		numberOfIterations()

		Returns the number of iterations performed by the algorithm.

		Returns
		-------
		int
			Number of iterations performed by the algorithm. Default: unlimited
		"""
		return (<_PageRank*>(self._this)).numberOfIterations()

	property norm:
		"""
		Property :code:`norm` can be one of the following:

		- networkit.centrality.Norm.l1norm
		- networkit.centrality.Norm.l2norm

		Set this property before calling run() to change norm computation.
		Default: networkit.centrality.Norm.l2norm
		"""
		def __get__(self):
			""" Get the norm used as stopping criterion. """
			return (<_PageRank*>(self._this)).norm
		def __set__(self, _Norm norm):
			""" Set the norm used as stopping criterion. """
			(<_PageRank*>(self._this)).norm = norm

	property maxIterations:
		"""
		Property :code:`maxIterations` sets a stopping criteria based on number
		of runs. Default: unlimited
		"""
		def __get__(self):
			""" Get the maximum number of iterations. """
			return (<_PageRank*>(self._this)).maxIterations
		def __set__(self, maxIterations):
			""" Set the maximum number of iterations. """
			if maxIterations < 0:
				raise Exception("Max iterations cannot be a negative number.")
			(<_PageRank*>(self._this)).maxIterations = maxIterations


cdef extern from "<networkit/centrality/SpanningEdgeCentrality.hpp>":

	cdef cppclass _SpanningEdgeCentrality "NetworKit::SpanningEdgeCentrality"(_Algorithm):
		_SpanningEdgeCentrality(_Graph G, double tol) except +
		void runApproximation() except +
		void runParallelApproximation() except +
		vector[double] scores() except +

cdef class SpanningEdgeCentrality(Algorithm):
	""" 
	SpanningEdgeCentrality(G, tol = 0.1)
	
	Computes the Spanning Edge centrality for the edges of the graph.

	Parameters
	----------
	G : networkit.Graph
		The graph.
	tol: float, optional
		Tolerance used for the approximation: with probability at least 1-1/n, the approximated 
		scores are within a factor 1+tol from the exact scores. Default: 0.1
	"""

	cdef Graph _G
	def __cinit__(self,  Graph G, double tol = 0.1):
		self._G = G
		self._this = new _SpanningEdgeCentrality(G._this, tol)

	def runApproximation(self):
		""" 
		runApproximation()

		Computes approximation of the Spanning Edge Centrality. This solves k linear systems, 
		where `k` is :math:`log(n)/(tol^2)`. The empirical running time is :math:`O(km)`, where n is the number of nodes
 		and m is the number of edges.
		"""
		return (<_SpanningEdgeCentrality*>(self._this)).runApproximation()

	def runParallelApproximation(self):
		""" 
		runParallelApproximation()
		
		Computes approximation (in parallel) of the Spanning Edge Centrality. This solves k linear systems, 
		where `k` is :math:`log(n)/(tol^2)`. The empirical running time is :math:`O(km)`, where n is the number of nodes
 	 	and m is the number of edges.
		"""
		return (<_SpanningEdgeCentrality*>(self._this)).runParallelApproximation()

	def scores(self):
		""" 
		scores()
		
		Get a vector containing the SEC score for each edge in the graph
		in ascending edge ID order.

		Returns
		-------
		list(float)
			The SEC scores.
		"""
		return (<_SpanningEdgeCentrality*>(self._this)).scores()

cdef extern from "<networkit/centrality/ApproxElectricalCloseness.hpp>":
	cdef cppclass _ApproxElectricalCloseness "NetworKit::ApproxElectricalCloseness"(_Centrality):
		_ApproxElectricalCloseness(_Graph G, double eps, double kappa) except +
		vector[double] getDiagonal() except +
		vector[double] computeExactDiagonal(double tol) except +

cdef class ApproxElectricalCloseness(Centrality):
	"""
	ApproxElectricalCloseness(G, eps = 0.1, kappa = 0.3)
	
	Approximates the electrical closeness of all the vertices of the graph by approximating the
	diagonal of the laplacian's pseudoinverse of G. Every element of the diagonal is
	guaranteed to have a maximum absolute error of eps. Based on "Approximation of the
	Diagonal of a Laplacian’s Pseudoinverse for Complex Network Analysis", Angriman et al., ESA
	2020. The algorithm does two steps: solves a linear system and samples uniform spanning trees
	(USTs). The parameter kappa balances the tolerance of solver for the linear system and the
	number of USTs to be sampled. A high value of kappa raises the tolerance (solver converges
	faster) but more USTs need to be sampled, vice versa for a low value of kappa.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	eps : float, optional
		Maximum absolute error of the elements in the diagonal.
	kappa : float, optional
		Balances the tolerance of the solver for the linear system and the
		number of USTs to be sampled.
	"""

	def __cinit__(self, Graph G, double eps = 0.1, double kappa = 0.3):
		self._G = G
		self._this = new _ApproxElectricalCloseness(G._this, eps, kappa)

	def getDiagonal(self):
		"""
		getDiagonal()

		Return an epsilon-approximation of the diagonal of the laplacian's pseudoinverse.

		Returns
		-------
		list(float)
			Approximation of the diagonal of the laplacian's pseudoinverse.
		"""
		return (<_ApproxElectricalCloseness*>self._this).getDiagonal()

	def computeExactDiagonal(self, double tol = 1e-9):
		"""
		computeExactDiagonal(tol = 1e-9)

		Compute and return the nearly-exact values of the diagonal of the laplacian's pseudoinverse.
		The values are computed by solving :math:`Lx = e_u - 1 / n` for every vertex `u` of the graph with a
		LAMG solver.

		Parameters
		----------
		tol : float
			Tolerance for the LAMG solver.

		Returns
		-------
		list(float)
			Nearly-exact values of the diagonal of the laplacian's pseudoinverse.
		"""
		return (<_ApproxElectricalCloseness*>self._this).computeExactDiagonal(tol)


cdef extern from "<networkit/centrality/ForestCentrality.hpp>":
	cdef cppclass _ForestCentrality "NetworKit::ForestCentrality"(_Centrality):
		_ForestCentrality(_Graph G, node root, double eps, double kappa) except +
		vector[double] &getDiagonal() except +

cdef class ForestCentrality(Centrality):
	"""
	ForestCentrality(G, root, eps=0.1, kappa=0.3)
	
	Approximates the forest closeness centrality of all the vertices of a graph by approximating
	the diagonal of the forest matrix of @a G. Every element of the diagonal is guaranteed to
	have a maximum absolute error of @a epsilon. Based on "New Approximation Algorithms for
	Forest Closeness Centrality - for Individual Vertices and Vertex Groups", van der Grinten et
	al, SDM 2021.
	The algorithm runs in two steps: (i) solving a linear system and (ii) sampling uniform
	spanning trees (USTs). The parameter @a kappa balances the tolerance of the linear sytem
	solver and the number of USTs to be sampled. A high value of @a kappa raises the tolerance
	(solver converges faster) but more USTs need to be sampled, vice versa for a low value of @a
	kappa. Note: the algorithm requires an augmented graph as input (see the reference paper for
	details). An augmented graphs can be generated with graphtools.createAugmentedGraph.

	Parameters
	----------
	G : networkit.Graph
		The input graph. Must be an augmented graph; an augmented graph can be crated with
		graphtools.createAugmentedGraph.
	root : int
		Root node of the augmented graph.
	epsilon : float, optional
		Maximum absolute error of the elements in the diagonal. Default: 0.1
	kappa : float, optional
		Balances the tolerance of the linear system solver and the number of USTs to be
		sampled. Default: 0.3
	"""

	def __cinit__(self, Graph G, node root, double eps = 0.1, double kappa = 0.3):
		self._G = G
		self._this = new _ForestCentrality(G._this, root, eps, kappa)

	def getDiagonal(self):
		"""
		getDiagonal()

		Return an epsilon-approximation of the diagonal of the forest matrix.

		Returns
		-------
		list(float)
			Approximation of the diagonal of the laplacian's pseudoinverse.
		"""
		return (<_ForestCentrality*>self._this).getDiagonal()

cdef extern from "<networkit/centrality/ApproxSpanningEdge.hpp>":
	cdef cppclass _ApproxSpanningEdge "NetworKit::ApproxSpanningEdge"(_Algorithm):
		_ApproxSpanningEdge(_Graph G, double eps) except +
		vector[edgeweight] scores() except +

cdef class ApproxSpanningEdge(Algorithm):
	"""
	ApproxSpanningEdge(G, eps = 0.1)

	Computes an epsilon-approximation of the spanning edge centrality of every edge of the input
	graph with probability :math:`(1 - 1/n)`, based on "Efficient Algorithms for Spanning Tree
	Centrality", Hayashi et al., IJCAI, 2016. This implementation also supports multi-threading.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	eps : float, optional
		Maximum additive error. Default: 0.1
	"""

	cdef Graph _G

	def __cinit__(self, Graph G, double eps = 0.1):
		self._G = G
		self._this = new _ApproxSpanningEdge(G._this, eps)

	def scores(self):
		"""
		scores()

		Return the spanning edge approximation for each edge of the graph
		in ascending edge ID order.

		Returns
		-------
		list
			Spanning edge approximation for each edge of the input graph.
		"""
		return (<_ApproxSpanningEdge*>(self._this)).scores();

def ranking(G, algorithm=Betweenness, normalized=False):
	""" 
	ranking(G, algorithm=networkit.centrality.Betweenness, normalized=False)

	Return a ranking of nodes by the specified centrality type
	
	Parameters
	----------
	G : networkit.Graph
		The input graph.
	algorithm : networkit.centrality.Centrality, optional
		Instance of centrality algorithm to run. Default: networkit.centrality.Betweenness
	normalized : bool, optional
		Set whether the ranking values should be normalized. Default: False

	Returns
	-------
	list(int)
		Ranking for nodes according to centrality algorithm.
	"""
	# FIXME: some centrality algorithms take more Parameters:
	centrality = algorithm(G, normalized)
	centrality.run()
	return centrality.ranking()

def scores(G, algorithm=Betweenness, normalized=False):
	""" 
	scores(G, algorithm=Betweenness, normalized=False)
	
	Return the centrality scores of nodes using the specified centrality type
	
	Parameters
	----------
	G : networkit.Graph
		The input graph.
	algorithm : networkit.centrality.Centrality, optional
		Instance of centrality algorithm to run. Default: networkit.centrality.Betweenness
	normalized : bool, optional
		Set whether the ranking values should be normalized. Default: False

	Returns
	-------
	list(int)
		Scores for nodes according to centrality algorithm.	
	"""
	centrality = algorithm(G, normalized)
	centrality.run()
	return centrality.scores()

def rankPerNode(ranking):
	"""
	rankPerNode(ranking)

	Returns ranks of all nodes sorted by their ID.

	Parameters
	----------
 	ranking: list(tuple(int, float))
		Ordered list of tuples (node, score).

	Returns
	-------
	list(int)
		For each node (sorted by node ID), the ranking of the node
	"""
	n_nodes = len(ranking)
	ranking_id = [0]*n_nodes
	for index, pair in enumerate(ranking):
		ranking_id[pair[0]] = index
	#we assign to all nodes the ranking of the first node with the same score
	for index, pair in enumerate(ranking):
			if index == 0:
				continue
			if pair[1] == ranking[index-1][1]:
				prev_node = ranking[index-1][0]
				ranking_id[pair[0]] = ranking_id[prev_node]
	return ranking_id

def relativeRankErrors(rx, ry):
	"""
	relativeRankErrors(rx, ry)

	Let :math:`r_x(u)` be the rank of node `u` in ranking `x`.
	The relative rank error of node `u` is defined as :math:`r_x(u) / r_y(u)`.

	Parameters
	----------
	rx : list(tuple(int, float))
		Ranking - ordered list of tuples (node, score).

	ry : list(tuple(int, float))
		Ranking - ordered list of tuples (node, score).

	Returns
	-------
	list(int)
		Rank errors ordered by node ID
	"""
	diff = []
	n = len(rx)
	if not(n == len(ry)):
		return diff
	rnode_x = rankPerNode(rx)
	rnode_y = rankPerNode(ry)
	for i in range(n):
		diff.append((rnode_x[i]+1)/(rnode_y[i]+1))
	return diff

class SpectralCentrality:
	"""
	SpectralCentrality(G, normalized=False)

	Abstract class to compute the spectral centrality of a graph. This class needs to be supplied with methods
	to generate the correct matrices and do the correct normalization.

	Parameters
	----------
	G : networkit.Graph
		The graph of which to compute the centrality.
	normalized : bool, optional
		Whether to normalize the results or not. Default: False
	"""
	def __init__(self, G, normalized=False):
		super(SpectralCentrality, self).__init__()

		self.graph = G
		self.normalized = normalized

		self.scoreList = None
		self.rankList = None
		self.evz = {}

	def prepareSpectrum(self):
		""" 
		prepareSpectrum()

		Not implemented yet.
		"""
		raise NotImplemented

	def normFactor(self):
		""" 
		normFactor()

		Not implemented yet.
		"""
		raise NotImplemented

	def run(self):
		""" 
		run()

		Runs computation of spectral centrality.
		"""
		self.prepareSpectrum()

		self.scoreList = None
		self.rankList = None
		self.evz = {}

		if self.normalized:
			normFactor = self.normFactor()
		else:
			normFactor = 1

		for v in self.graph.iterNodes():
			self.evz[v] = self.eigenvector[v] * normFactor
		return self

	def scores(self):
		""" 
		scores()
		
		Get a vector containing the SpectralCentrality score.

		Returns
		-------
		list(float)
			The SpectralCentrality scores.
		"""
		if self.scoreList is None:
			self.scoreList = [v for k,v in self.evz.items()]

		return self.scoreList

	def ranking(self):
		""" 
		ranking()

		Return a ranking of nodes by SpectralCentrality.

		Returns
		-------
		list(float)
			Ranking for nodes according to SpectralCentrality.
		"""
		if self.rankList is None:
			self.rankList = sorted(self.evz.items(),key=lambda x: float(x[1]), reverse=True)
		return self.rankList


class SciPyEVZ(SpectralCentrality):
	"""
	SciPyEVZ(G, normalized=False)

	Compute Eigenvector centrality using algebraic meh

	Parameters
	----------
	G : networkit.Graph
		The graph of which to compute the centrality.
	normalized : bool, optional
		Whether to normalize the results or not. Default: False
	"""
	def __init__(self, G, normalized=False):
		if G.isDirected():
			raise NotImplementedError("Not implemented for directed graphs; use centrality.EigenvectorCentrality instead")
		super(SciPyEVZ, self).__init__(G, normalized=normalized)

	def _length(self, vector):
		square = sum([val * val for val in vector])
		return math.sqrt(square)

	def normFactor(self):
		"""
		normFactor()

		Returns the norm factor.

		Returns
		-------
		float
			Norm factor.
		"""
		return 1 / self._length(self.eigenvector)

	def prepareSpectrum(self):
		"""
		prepareSpectrum()

		Prepare the computation of SciPyEVZ.
		"""
		spectrum = adjacencyEigenvector(self.graph, i=0)
		self.eigenvector = spectrum[1]
		self.eigenvalue = spectrum[0]

class SciPyPageRank(SpectralCentrality):
	"""
	SciPyPageRank(G, damp=0.95, normalized=False)

	Compute Eigenvector centrality using algebraic meh

	Parameters
	----------
	G : networkit.Graph
		The graph of which to compute the centrality.
	damp : float, optional
		Damping factor used for computation. Default: 0.95
	normalized : bool, optional
		Whether to normalize the results or not. Default: False
	"""
	def __init__(self, G, damp=0.95, normalized=False):
		super(SciPyPageRank, self).__init__(G, normalized=normalized)

		self.damp = damp

	def _length(self, vector):
		return sum(vector)

	def normFactor(self):
		"""
		normFactor()

		Returns the norm factor.

		Returns
		-------
		float
			Norm factor.
		"""
		return 1 / self._length(self.eigenvector)

	def prepareSpectrum(self):
		"""
		prepareSpectrum()

		Prepare the computation of SciPyPageRank.
		"""
		prMatrix = PageRankMatrix(self.graph, self.damp)
		spectrum = symmetricEigenvectors(prMatrix, cutoff=0, reverse=False)

		self.eigenvector = spectrum[1][0]
		self.eigenvalue = spectrum[0][0]
