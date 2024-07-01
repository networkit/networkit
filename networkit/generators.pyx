# distutils: language=c++

import math
import subprocess
import os
import tempfile
import scipy

from libcpp.vector cimport vector
from libcpp cimport bool as bool_t
from libcpp.utility cimport pair

from .base cimport _Algorithm, Algorithm
from .centrality import DegreeCentrality, LocalPartitionCoverage
from .community import PLM
from .dynamics cimport _GraphEvent, GraphEvent
from .graph cimport _Graph, Graph
from .graphtools import GraphTools
from .structures cimport _Partition, Partition, count, index, node, coordinate

cdef extern from "<networkit/viz/Point.hpp>" namespace "NetworKit" nogil:

	cdef cppclass Point[T]:
		Point()
		Point(T x, T y)
		T& operator[](const index i) except +
		T& at(const index i) except +

	cdef cppclass _Point2D "NetworKit::Point2D":
		_Point2D()
		pair[coordinate, coordinate] asPair()

cdef object toPoint2DVector(const vector[_Point2D]& v):
	return [v[i].asPair() for i in range(v.size())]

cdef object toNodePoint2DVector(const vector[pair[node, _Point2D]]& v):
	return [(v[i].first, v[i].second.asPair()) for i in range(v.size())]

cdef extern from "<networkit/generators/StaticGraphGeneratorBase.hpp>":

	cdef cppclass _StaticGraphGenerator "NetworKit::StaticGraphGenerator":
		_StaticGraphGenerator()
		_Graph generate() except +

cdef class StaticGraphGenerator:
	""" Abstract base class for static graph generators """
	cdef _StaticGraphGenerator *_this

	def __init__(self, *args, **namedargs):
		if type(self) == StaticGraphGenerator:
			raise RuntimeError("Error, you may not use StaticGraphGenerator directly, use a sub-class instead")

	def __cinit__(self, *args, **namedargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

	def generate(self):
		"""
		generate()
		
		Generates the graph.

		Returns
		-------
		networkit.Graph
			The generated graph.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return Graph().setThis(self._this.generate())

cdef extern from "<networkit/generators/BarabasiAlbertGenerator.hpp>":

	cdef cppclass _BarabasiAlbertGenerator "NetworKit::BarabasiAlbertGenerator"(_StaticGraphGenerator):
		_BarabasiAlbertGenerator() except +
		_BarabasiAlbertGenerator(count k, count nMax, count n0, bool_t batagelj) except +
		_BarabasiAlbertGenerator(count k, count nMax, const _Graph & initGraph, bool_t batagelj) except +

cdef class BarabasiAlbertGenerator(StaticGraphGenerator):
	"""
	BarabasiAlbertGenerator(k, nMax, n0=0, batagelj=True)

	This generator implements the preferential attachment model as introduced by Barabasi and Albert[1].
	The original algorithm is very slow and thus, the much faster method from Batagelj and Brandes[2] is
	implemented and the current default.
	The original method can be chosen by setting \p batagelj to false.
	[1] Barabasi, Albert: Emergence of Scaling in Random Networks http://arxiv.org/pdf/cond-mat/9910332.pdf
	[2] ALG 5 of Batagelj, Brandes: Efficient Generation of Large Random Networks https://kops.uni-konstanz.de/bitstream/handle/123456789/5799/random.pdf?sequence=1

	Parameters
	----------
	k : int
		Number of edges that come with a new node.
	nMax : int
		Maximum number of nodes produced.
	n0 : int or networkit.Graph
		Number of starting nodes or the initial starting graph. Default: 0
	batagelj : bool
		Specifies whether to use batagelj's method or the original one. Default: True
	"""

	def __cinit__(self, count k, count nMax, n0=0, bool_t batagelj=True):
		if isinstance(n0, Graph):
			self._this = new _BarabasiAlbertGenerator(k, nMax, (<Graph>n0)._this, batagelj)
		else:
			self._this = new _BarabasiAlbertGenerator(k, nMax, <count>n0, batagelj)

	@classmethod
	def fit(cls, Graph G, scale=1):
		"""
		fit(G, scale=1)

		Fit model to input graph

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		scale : int, optional
			Scale the maximum number of nodes by a factor. Default: 1

		Returns
		-------
		networkit.Graph
			A new scaled graph.
		"""
		(n, m) = GraphTools.size(G)
		k = math.floor(m / n)
		return cls(nMax=scale * n, k=k, n0=k)

cdef extern from "<networkit/generators/PubWebGenerator.hpp>":

	cdef cppclass _PubWebGenerator "NetworKit::PubWebGenerator"(_StaticGraphGenerator):
		_PubWebGenerator(count numNodes, count numberOfDenseAreas, float neighborhoodRadius, count maxNumberOfNeighbors) except +
		const vector[_Point2D]& getCoordinates()


cdef class PubWebGenerator(StaticGraphGenerator):
	""" 
	PubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)
	
	Generates a static graph that resembles an assumed geometric distribution of nodes in
	a P2P network.

	The basic structure is to distribute points randomly in the unit torus
	and to connect vertices close to each other (at most @a neighRad distance and none of
	them already has @a maxNeigh neighbors). The distribution is chosen to get some areas with
	high density and others with low density. There are @a numDenseAreas dense areas, which can
	overlap. Each area is circular, has a certain position and radius and number of points.
	These values are strored in @a denseAreaXYR and @a numPerArea, respectively.

	Used and described in more detail in J. Gehweiler, H. Meyerhenke: A Distributed
	Diffusive Heuristic for Clustering a Virtual P2P Supercomputer. In Proc. 7th High-Performance
	Grid Computing Workshop (HPGC'10), in conjunction with 24th IEEE Internatl. Parallel and
	Distributed Processing Symposium (IPDPS'10), IEEE, 2010.

	Parameters
	----------
	numNodes : int
		Up to a few thousand (possibly more if visualization is not desired and quadratic
		time complexity has been resolved)
	numberOfDenseAreas : int
		Depending on number of nodes, e.g. [8, 50]
	neighborhoodRadius : float
		The higher, the better the connectivity [0.1, 0.35]
	maxNumberOfNeighbors : int
		Maximum degree, a higher value corresponds to better connectivity [4, 40]
	"""

	def __cinit__(self, numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors):
		self._this = new _PubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	def getCoordinates(self):
		"""
		getCoordinates()
		
		Returns a list of coordinates
		
		Returns
		-------
		list(tuple(float,float))
			2D coordinates of all nodes in the graph.
		"""
		return toPoint2DVector((<_PubWebGenerator*>(self._this)).getCoordinates())

cdef extern from "<networkit/generators/DynamicPubWebGenerator.hpp>":
	cdef cppclass _DynamicPubWebGenerator "NetworKit::DynamicPubWebGenerator":
		_DynamicPubWebGenerator(count numNodes, count numberOfDenseAreas,
			float neighborhoodRadius, count maxNumberOfNeighbors) except +
		vector[_GraphEvent] generate(count nSteps) except +
		_Graph getGraph() except +
		vector[_Point2D] getCoordinates()
		vector[pair[node, _Point2D]] getNewCoordinates()

cdef class DynamicPubWebGenerator:
	"""
	DynamicPubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	Dynamic variant of networkit.generators.PubWebGenerator.

	Parameters
	----------
	numNodes : int
		Up to a few thousand (possibly more if visualization is not desired and quadratic
		time complexity has been resolved)
	numberOfDenseAreas : int
		Depending on number of nodes, e.g. [8, 50]
	neighborhoodRadius : float
		The higher, the better the connectivity [0.1, 0.35]
	maxNumberOfNeighbors : int
		Maximum degree, a higher value corresponds to better connectivity [4, 40]
	"""
	cdef _DynamicPubWebGenerator* _this

	def __cinit__(self, numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors):
		self._this = new _DynamicPubWebGenerator(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors)

	def __dealloc__(self):
		del self._this

	def generate(self, nSteps):
		""" 
		generate(nSteps)
		
		Generate event stream.

		Parameters
		----------
		nSteps : int
			Number of time steps in the event stream.

		Returns
		-------
		list(networkit.dynamics.GraphEvent)
			List of graph events.
		"""
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]

	def getGraph(self):
		"""
		getGraph()

		Returns current graph

		Returns
		-------
		networkit.Graph
			The resulting graph.
		"""
		return Graph().setThis(self._this.getGraph())

	def getCoordinates(self):
		"""
		getCoordinates()
		
		Returns a list of coordinates from the current state.
		
		Returns
		-------
		list(tuple(float,float))
			2D coordinates of all nodes in the graph.
		"""
		return toPoint2DVector((<_DynamicPubWebGenerator*>(self._this)).getCoordinates())

	def getNewCoordinates(self):
		"""
		getNewCoordinates()
		
		Get list of nodes and coordinates of points added during last generate call.
		
		Returns
		-------
		list(int, tuple(float,float))
			List of node ids and corresponding coordinates.
		"""
		return toNodePoint2DVector((<_DynamicPubWebGenerator*>(self._this)).getNewCoordinates())

cdef extern from "<networkit/generators/ErdosRenyiGenerator.hpp>":

	cdef cppclass _ErdosRenyiGenerator "NetworKit::ErdosRenyiGenerator"(_StaticGraphGenerator):
		_ErdosRenyiGenerator(count nNodes, double prob, bool_t directed, bool_t selfLoops) except +

cdef class ErdosRenyiGenerator(StaticGraphGenerator):
	"""
	ErdosRenyiGenerator(count nNodes, double prob, directed = False, selfLoops = False)
	
	Creates random graphs in the G(n,p) model.
	The generation follows Vladimir Batagelj and Ulrik Brandes: "Efficient
	generation of large random networks", Phys Rev E 71, 036113 (2005).
	
	Creates G(nNodes, prob) graphs.

	Parameters
	----------
	nNodes : int
		Number of nodes n in the graph.
	prob : float
		Probability of existence for each edge p.
	directed : bool, optional
		Generates a directed graph. Default: False.
	selfLoops : bool, optional
		Allows self-loops to be generated (only for directed graphs). Default: False.
	"""

	def __cinit__(self, nNodes, prob, directed = False, selfLoops = False):
		self._this = new _ErdosRenyiGenerator(nNodes, prob, directed, selfLoops)

	@classmethod
	def fit(cls, Graph G, scale=1):
		"""
		fit(G, scale=1)

		Fit model to input graph

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		scale : int, optional
			Scale the maximum number of nodes by a factor. Default: 1

		Returns
		-------
		networkit.Graph
			A new scaled graph.
		"""
		(n, m) = GraphTools.size(G)
		if G.isDirected():
			raise Exception("TODO: figure out scaling scheme for directed graphs")
		else:
			p = (2 * m) / (scale * n * (n-1))
		return cls(scale * n, p)

cdef extern from "<networkit/generators/DorogovtsevMendesGenerator.hpp>":

	cdef cppclass _DorogovtsevMendesGenerator "NetworKit::DorogovtsevMendesGenerator"(_StaticGraphGenerator):
		_DorogovtsevMendesGenerator(count nNodes) except +

cdef class DorogovtsevMendesGenerator(StaticGraphGenerator):
	""" 
	DorogovtsevMendesGenerator(nNodes)

	Generates a graph according to the Dorogovtsev-Mendes model.

	Parameters
	----------
	nNodes : int
		Number of nodes in the target graph.
	"""

	def __cinit__(self, nNodes):
		self._this = new _DorogovtsevMendesGenerator(nNodes)

	@classmethod
	def fit(cls, Graph G, scale=1):
		"""
		fit(G, scale=1)

		Fit model to input graph

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		scale : int, optional
			Scale the maximum number of nodes by a factor. Default: 1

		Returns
		-------
		networkit.Graph
			A new scaled graph.
		"""
		return cls(scale * G.numberOfNodes())

cdef extern from "<networkit/generators/ClusteredRandomGraphGenerator.hpp>":

	cdef cppclass _ClusteredRandomGraphGenerator "NetworKit::ClusteredRandomGraphGenerator"(_StaticGraphGenerator):
		_ClusteredRandomGraphGenerator(count, count, double, double) except +
		_Partition getCommunities() except +

cdef class ClusteredRandomGraphGenerator(StaticGraphGenerator):
	""" 
	ClusteredRandomGraphGenerator(count, count, pin, pout)
	
	The ClusteredRandomGraphGenerator class is used to create a clustered random graph.

	The number of nodes and the number of edges are adjustable as well as the probabilities
	for intra-cluster and inter-cluster edges.

	In parallel the generated graph is not deterministic. To ensure determinism, use a single thread.

	Parameters
	----------
	n : int
		Number of nodes.
	k : int
		Number of clusters.
	pin : float
		Intra-cluster edge probability.
	pout : float
		Inter-cluster edge probability.
	"""

	def __cinit__(self, n, k, pin, pout):
		self._this = new _ClusteredRandomGraphGenerator(n, k, pin, pout)

	def getCommunities(self):
		""" 
		getCommunities()
		
		Returns the generated ground truth clustering.

		Returns
		-------
		networkit.Partition
			The generated ground truth clustering.
		"""
		return Partition().setThis((<_ClusteredRandomGraphGenerator*>(self._this)).getCommunities())

cdef extern from "<networkit/generators/ChungLuGenerator.hpp>":

	cdef cppclass _ChungLuGenerator "NetworKit::ChungLuGenerator"(_StaticGraphGenerator):
		_ChungLuGenerator(vector[count] degreeSequence) except +

cdef class ChungLuGenerator(StaticGraphGenerator):
	"""
	ChungLuGenerator(degreeSequence)

	Given an arbitrary degree sequence, the Chung-Lu generative model
	will produce a random graph with the same expected degree sequence.

	see Chung, Lu: The average distances in random graphs with given expected degrees
	and Chung, Lu: Connected Components in Random Graphs with Given Expected Degree Sequences.
	Aiello, Chung, Lu: A Random Graph Model for Massive Graphs describes a different generative model
	which is basically asymptotically equivalent but produces multi-graphs.

	Parameters
	----------
	degreeSequence : list(float)
		Input degree sequence used to generate the graph.
	"""

	def __cinit__(self, vector[count] degreeSequence):
		self._this = new _ChungLuGenerator(degreeSequence)

	@classmethod
	def fit(cls, Graph G, scale=1):
		"""
		fit(G, scale=1)

		Fit model to input graph

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		scale : int, optional
			Scale the maximum number of nodes by a factor. Default: 1

		Returns
		-------
		networkit.Graph
			A new scaled graph.
		"""
		(n, m) = GraphTools.size(G)
		degSeq = DegreeCentrality(G).run().scores()
		return cls(degSeq * scale)

cdef extern from "<networkit/generators/HyperbolicGenerator.hpp>":

	cdef cppclass _HyperbolicGenerator "NetworKit::HyperbolicGenerator"(_StaticGraphGenerator):
		_HyperbolicGenerator(count nodes,  double k, double gamma, double T) except +
		void setLeafCapacity(count capacity) except +
		void setTheoreticalSplit(bool_t split) except +
		void setBalance(double balance) except +
		vector[double] getElapsedMilliseconds() except +
		_Graph generate(vector[double] angles, vector[double] radii, double R, double T) except +

cdef class HyperbolicGenerator(StaticGraphGenerator):
	"""
	HyperbolicGenerator(n, k=6, gamma=3, T=0)
	
	The Hyperbolic Generator distributes points in hyperbolic space and adds edges between points 
	with a probability depending on their distance. The resulting graphs have a power-law degree 
	distribution, small diameter and high clustering coefficient.
	For a temperature of 0, the model resembles a unit-disk model in hyperbolic space.

	Parameters
	----------
	n : int
		Number of nodes.
	k : float, optional
		Average degree. Default: 6.0
	gamma : float, optional
		Exponent of power-law degree distribution. Default: 3.0
	T : float, optional
		Temperature of statistical model. Default: 0.0
	"""

	def __cinit__(self,  n, k=6, gamma=3, T=0):
		if gamma <= 2:
				raise ValueError("Exponent of power-law degree distribution must be > 2")
		self._this = new _HyperbolicGenerator(n, k, gamma, T)

	def setLeafCapacity(self, capacity):
		"""
		setLeafCapacity(capacity)
		
		Set the capacity of a quadtree leaf.

		Parameters
		----------
		capacity : int
			Tuning parameter.
		"""
		(<_HyperbolicGenerator*>(self._this)).setLeafCapacity(capacity)

	def setBalance(self, balance):
		"""
		setBalance(balance)

		Set the balance of the quadtree. Value should be between 0.0 and 1.0.
		Without modification this is set to 0.5.

		Parameters
		----------
		balance : float
			Balance factor between 0.0 and 1.0.
		"""
		(<_HyperbolicGenerator*>(self._this)).setBalance(balance)

	def setTheoreticalSplit(self, theoreticalSplit):
		"""
		setTheoreticalSplit(theoreticalSplit)

		When using a theoretically optimal split, the quadtree will be flatter, 
		but running time usually longer.

		Parameters
		----------
		theoreticalSplit : bool
			Whether to use the theoretically optimal split. Default: False.
		"""
		(<_HyperbolicGenerator*>(self._this)).setTheoreticalSplit(theoreticalSplit)

	def getElapsedMilliseconds(self):
		"""
		getElapsedMilliseconds()

		Get running time of generator for each.

		Returns
		-------
		list(float)
			Running time of the generator for each thread.
		"""
		return (<_HyperbolicGenerator*>(self._this)).getElapsedMilliseconds()

	def generate_advanced(self, angles, radii, R, T=0):
		"""
		generate_advanced(angles, radii, R, T=0)

		Generate a graph with overriding parameters.

		Parameters
		----------
		angles : list(float)
			List with angles of node positions.
		radii : list(float)
			List with radii of node positions.
		R : float
			Radius of poincare disk to place nodes in.
		T : float, optional
			Edges are added for nodes closer to each other than threshold T. Default: 0.0
		"""
		return Graph(0).setThis((<_HyperbolicGenerator*>(self._this)).generate(angles, radii, R, T))

	@classmethod
	def fit(cls, Graph G, scale=1):
		"""
		fit(G, scale=1)

		Fit model to input graph

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		scale : int, optional
			Scale the maximum number of nodes by a factor. Default: 1

		Returns
		-------
		networkit.Graph
			A new scaled graph.
		"""
		degSeq = DegreeCentrality(G).run().scores()
		gamma = max(-1 * PowerlawDegreeSequence(degSeq).getGamma(), 2.1)
		(n, m) = GraphTools.size(G)
		k = 2 * (m / n)
		return cls(n * scale, k, gamma)

cdef extern from "<networkit/generators/PowerlawDegreeSequence.hpp>":

	cdef cppclass _PowerlawDegreeSequence "NetworKit::PowerlawDegreeSequence":
		_PowerlawDegreeSequence(count minDeg, count maxDeg, double gamma) except +
		_PowerlawDegreeSequence(_Graph) except +
		_PowerlawDegreeSequence(vector[double]) except +
		void setMinimumFromAverageDegree(double avgDeg) nogil except +
		void setGammaFromAverageDegree(double avgDeg, double minGamma, double maxGamma) nogil except +
		double getExpectedAverageDegree() except +
		count getMinimumDegree() const
		count getMaximumDegree() const
		double getGamma() const
		double setGamma(double) const
		void run() nogil except +
		vector[count] getDegreeSequence(count numNodes) except +
		count getDegree() except +

cdef class PowerlawDegreeSequence:
	"""
	PowerlawDegreeSequence(minDeg, maxDeg, gamma)

	Other calling possibilities:

	- PowerlawDegreeSequence(G)
	- PowerlawDegreeSequence(degreeSequence)

	Generates a powerlaw degree sequence with the given minimum and maximum degree, the powerlaw exponent gamma

	If a list of degrees or a graph is given instead of a minimum degree, the class uses the minimum and maximum
	value of the sequence and fits the exponent such that the expected average degree is the actual average degree.

	Parameters
	----------
	minDeg : int
		The minium degree.
	maxDeg : int, optional
		The maximum degree.
	gamma : float, optional
		The powerlaw exponent. Default: -2.0
	G : networkit.Graph, alternative
		The input graph.
	degreeSequence : list(int), alternative
		List of degrees to fit.
	"""
	cdef _PowerlawDegreeSequence *_this

	def __cinit__(self, minDeg, count maxDeg = 0, double gamma = -2):
		if isinstance(minDeg, Graph):
			self._this = new _PowerlawDegreeSequence((<Graph>minDeg)._this)
		else:
			try:
				self._this = new _PowerlawDegreeSequence(<vector[double]?>minDeg)
			except TypeError:
				self._this = new _PowerlawDegreeSequence((<count?>minDeg), maxDeg, gamma)

	def __dealloc__(self):
		del self._this

	def setMinimumFromAverageDegree(self, double avgDeg):
		"""
		setMinimumFromAverageDegree(avgDeg)
		
		Tries to set the minimum degree such that the specified average degree is expected.

		Parameters
		----------
		avgDeg : float
			The average degree that shall be approximated.
		"""
		with nogil:
			self._this.setMinimumFromAverageDegree(avgDeg)
		return self

	def setGammaFromAverageDegree(self, double avgDeg, double minGamma = -1, double maxGamma = -6):
		"""
		setGammaFromAverageDegree(avgDeg, minGamma = -1, maxGamma = -6)

		Tries to set the powerlaw exponent gamma such that the specified average degree is expected.

		Parameters
		----------
		avgDeg : float
			The average degree that shall be approximated.
		minGamma : float, optional
			The minimum gamma to use. Default: -1.0
		maxGamma : float, optional
			The maximum gamma to use. Default: -6.0
		"""
		with nogil:
			self._this.setGammaFromAverageDegree(avgDeg, minGamma, maxGamma)
		return self

	def getExpectedAverageDegree(self):
		"""
		getExpectedAverageDegree()

		Returns the expected average degree. Note: run needs to be called first.

		Returns
		-------
		float
			The expected average degree.
		"""
		return self._this.getExpectedAverageDegree()

	def getMinimumDegree(self):
		"""
		getMinimumDegree()

		Returns the minimum degree.

		Returns
		-------
		int
			The minimum degree.
		"""
		return self._this.getMinimumDegree()

	def setGamma(self, double gamma):
		"""
		setGamma(gamma)

		Set the exponent gamma.

		Parameters
		----------
		gamma : float
			The exponent to set.
		"""
		self._this.setGamma(gamma)
		return self

	def getGamma(self):
		"""
		getGamma()

		Get the exponent gamma.

		Returns
		-------
		float
			The exponent gamma.
		"""
		return self._this.getGamma()

	def getMaximumDegree(self):
		"""
		getMaximumDegree()

		Get the maximum degree

		Returns
		-------
		int
			The maximum degree
		"""
		return self._this.getMaximumDegree()

	def run(self):
		"""
		run()

		Executes the generation of the probability distribution.
		"""
		with nogil:
			self._this.run()
		return self

	def getDegreeSequence(self, count numNodes):
		"""
		getDegreeSequence(numNodes)

		Returns a degree sequence with even degree sum.

		Parameters
		----------
		numNodes : int
			The number of nodes/degrees that shall be returned.

		Returns
		-------
		list(int)
			The generated degree sequence.
		"""
		return self._this.getDegreeSequence(numNodes)

	def getDegree(self):
		"""
		getDegree()

		Returns a degree drawn at random with a power law distribution.

		Returns
		-------
		int
			The generated random degree.
		"""
		return self._this.getDegree()

cdef extern from "<networkit/generators/LFRGenerator.hpp>":

	cdef cppclass _LFRGenerator "NetworKit::LFRGenerator"(_Algorithm):
		_LFRGenerator(count n) except +
		void setDegreeSequence(vector[count] degreeSequence) nogil except +
		void generatePowerlawDegreeSequence(count avgDegree, count maxDegree, double nodeDegreeExp) nogil except +
		void setCommunitySizeSequence(vector[count] communitySizeSequence) nogil except +
		void setPartition(_Partition zeta) nogil except +
		void generatePowerlawCommunitySizeSequence(count minCommunitySize, count maxCommunitySize, double communitySizeExp) nogil except +
		void setMu(double mu) nogil except +
		void setMu(const vector[double] & mu) nogil except +
		void setMuWithBinomialDistribution(double mu) nogil except +
		_Graph getGraph() except +
		_Partition getPartition() except +
		_Graph generate() except +

cdef class LFRGenerator(Algorithm):
	"""
	LFRGenerator(n)

	The LFR clustered graph generator as introduced by Andrea Lancichinetti, Santo Fortunato, and Filippo Radicchi.

	The community assignment follows the algorithm described in
	"Benchmark graphs for testing community detection algorithms". The edge generation is however taken from their follow-up publication
	"Benchmarks for testing community detection algorithms on directed and weighted graphs with overlapping communities". Parts of the
	implementation follow the choices made in their implementation which is available at https://sites.google.com/site/andrealancichinetti/software
	but other parts differ, for example some more checks for the realizability of the community and degree size distributions are done
	instead of heavily modifying the distributions.

	The edge-switching markov-chain algorithm implementation in NetworKit is used which is different from the implementation in the original LFR benchmark.

	You need to set a degree sequence, a community size sequence and a mu using the additionally provided set- or generate-methods.

	Parameters
	----------
	n : int
		The number of nodes.
	"""
	params = {}
	paths = {}

	def __cinit__(self, count n):
		self._this = new _LFRGenerator(n)

	def setDegreeSequence(self, vector[count] degreeSequence):
		"""
		setDegreeSequence(degreeSequence)

		Set the given degree sequence.

		Parameters
		----------
		degreeSequence : list(int)
			The degree sequence that shall be used by the generator.
		"""
		with nogil:
			(<_LFRGenerator*>(self._this)).setDegreeSequence(degreeSequence)
		return self

	def generatePowerlawDegreeSequence(self, count avgDegree, count maxDegree, double nodeDegreeExp):
		"""
		generatePowerlawDegreeSequence(avgDegree, maxDegree, nodeDegreeExp)
		
		Generate and set a power law degree sequence using the given average and maximum degree with the given exponent.

		Parameters
		----------
		avgDegree : int
			The average degree of the created graph.
		maxDegree : int
			The maximum degree of the created graph.
		nodeDegreeExp : float
			The (negative) exponent of the power law degree distribution of the node degrees.
		"""
		with nogil:
			(<_LFRGenerator*>(self._this)).generatePowerlawDegreeSequence(avgDegree, maxDegree, nodeDegreeExp)
		return self

	def setCommunitySizeSequence(self, vector[count] communitySizeSequence):
		"""
		setCommunitySizeSequence(communitySizeSequence)

		Set the given community size sequence.

		Parameters
		----------
		communitySizeSequence : list(float)
			The community sizes that shall be used.
		"""
		with nogil:
			(<_LFRGenerator*>(self._this)).setCommunitySizeSequence(communitySizeSequence)
		return self

	def setPartition(self, Partition zeta not None):
		"""
		setPartition(zeta)

		Set the partition, this replaces the community size sequence and the random assignment of the nodes to communities.

		Parameters
		----------
		zeta : networkit.Partition
			The partition to use.
		"""
		with nogil:
			(<_LFRGenerator*>(self._this)).setPartition(zeta._this)
		return self

	def generatePowerlawCommunitySizeSequence(self, count minCommunitySize, count maxCommunitySize, double communitySizeExp):
		"""
		generatePowerlawCommunitySizeSequence(minCommunitySize, maxCommunitySize, communitySizeExp)
		
		Generate a powerlaw community size sequence with the given minimum and maximum size and the given exponent.

		Parameters
		----------
		minCommunitySize : int
			The minimum community size.
		maxCommunitySize : int
			The maximum community size.
		communitySizeExp : float
			The (negative) community size exponent of the power law degree distribution of the community sizes.
		"""
		with nogil:
			(<_LFRGenerator*>(self._this)).generatePowerlawCommunitySizeSequence(minCommunitySize, maxCommunitySize, communitySizeExp)
		return self

	def setMu(self, mu):
		"""
		setMu(mu)

		Set the mixing parameter, this is the fraction of neighbors of each node that do not belong to the node's own community.

		This can either be one value for all nodes or an iterable of values for each node.

		Parameters
		----------
		mu : float or list(float)
			The mixing coefficient(s), i.e. the factor of the degree that shall be inter-cluster degree
		"""
		try:
			(<_LFRGenerator*>(self._this)).setMu(<vector[double]?>mu)
		except TypeError:
			(<_LFRGenerator*>(self._this)).setMu(<double?>mu)
		return self

	def setMuWithBinomialDistribution(self, double mu):
		"""
		setMuWithBinomialDistribution(mu)

		Set the internal degree of each node using a binomial distribution such that the expected mixing parameter is the given @a mu.

		The mixing parameter is for each node the fraction of neighbors that do not belong to the node's own community.

		Parameters
		----------
		mu : float
			The expected mu that shall be used.
		"""
		with nogil:
			(<_LFRGenerator*>(self._this)).setMuWithBinomialDistribution(mu)
		return self

	def getGraph(self):
		"""
		getGraph()

		Return the generated Graph.

		Returns
		-------
		networkit.Graph
			The generated graph.
		"""
		return Graph().setThis((<_LFRGenerator*>(self._this)).getGraph())

	def generate(self, useReferenceImplementation=False):
		"""
		generate(useReferenceImplementation=False)

		Generates and returns the graph. Wrapper for the StaticGraphGenerator interface.

		Parameters
		----------
		useReferenceImplementation : bool
			Sets whether the reference implmentation should be used for generating. Default: False

		Returns
		-------
		networkit.Graph
			The generated graph.
		"""
		if useReferenceImplementation:
			from networkit import graphio
			os.system("{0}/benchmark {1}".format(self.paths["refImplDir"], self.params["refImplParams"]))
			return graphio.readGraph("network.dat", graphio.Format.EdgeListTabOne)
		return Graph().setThis((<_LFRGenerator*>(self._this)).generate())

	def getPartition(self):
		"""
		getPartition()

		Return the generated Partiton.

		Returns
		-------
		networkit.Partition
			The generated partition.
		"""
		return Partition().setThis((<_LFRGenerator*>(self._this)).getPartition())

	@classmethod
	def setPathToReferenceImplementationDir(cls, path):
		"""
		setPaths(path)

		Helper function for fit. Sets an output folder for writing the result to.

		Parameters
		----------
		path : str
			String containing a path.
		"""
		cls.paths["refImplDir"] = path


	@classmethod
	def fit(cls, Graph G, scale=1, vanilla=False, communityDetectionAlgorithm=PLM, plfit=False):
		""" 
		fit(G, scale=1, vanilla=False, communityDetectionAlgorithm=PLM, plfit=False)

		Fit model to input graph

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		scale : int, optional
			Scale the maximum number of nodes by a factor. Default: 1
		vanilla : bool, optional
			If set to True, fit power law to degree distribution. Otherwise fit to community sequence.
		communityDetectionAlgorithm : nk.community.CommunityDetector, optional
			Community detection algorithm used for fitting. Default: nk.community.PLM
		plfit : bool, optional
			If set to True, power law fitting is enabled. Default: False
		Returns
		-------
		networkit.Graph
			A new scaled graph.
		"""
		(n, m) = GraphTools.size(G)
		# detect communities
		communities = communityDetectionAlgorithm(G).run().getPartition()
		# get degree sequence
		degSeq = DegreeCentrality(G).run().scores()
		# set number of nodes
		gen = cls(n * scale)
		if vanilla:
			# fit power law to degree distribution and generate degree sequence accordingly
			#print("fit power law to degree distribution and generate degree sequence accordingly")
			avgDegree = int(sum(degSeq) / len(degSeq))
			maxDegree = max(degSeq)
			if plfit:
				degSeqGen = PowerlawDegreeSequence(G)
				nodeDegreeExp = -1 * degSeqGen.getGamma()
				degSeqGen.run()
				gen.setDegreeSequence(degSeqGen.getDegreeSequence(n * scale))
			else:
				nodeDegreeExp = 2
				gen.generatePowerlawDegreeSequence(avgDegree, maxDegree, -1 * nodeDegreeExp)
			print(avgDegree, maxDegree, nodeDegreeExp)
			# fit power law to community size sequence and generate accordingly
			#print("fit power law to community size sequence and generate accordingly")
			communitySize = communities.subsetSizes()
			communityAvgSize = int(sum(communitySize) / len(communitySize))
			communityMaxSize = max(communitySize)
			communityMinSize = min(communitySize)

			localCoverage = LocalPartitionCoverage(G, communities).run().scores()
			mu = 1.0 - sum(localCoverage) / len(localCoverage)
			# check if largest possible internal degree can fit in the largest possible community
			if math.ceil((1.0 - mu) * maxDegree) >= communityMaxSize:
				# Make the maximum community size 5% larger to make it more likely
				# the largest generated degree will actually fit.
				communityMaxSize = math.ceil(((1.0 - mu) * maxDegree + 1) * 1.05)
				print("Increasing maximum community size to fit the largest degree")

			if plfit:
				communityExp = -1 * PowerlawDegreeSequence(communityMinSize, communityMaxSize, -1).setGammaFromAverageDegree(communityAvgSize).getGamma()
			else:
				communityExp = 1
			pl = PowerlawDegreeSequence(communityMinSize, communityMaxSize, -1 * communityExp)

			try: # it can be that the exponent is -1 because the average would be too low otherwise, increase minimum to ensure average fits.
				pl.setMinimumFromAverageDegree(communityAvgSize)
				communityMinSize = pl.getMinimumDegree()
			except RuntimeError: # if average is too low with chosen exponent, this might not work...
				pl.run()
				print("Could not set desired average community size {}, average will be {} instead".format(communityAvgSize, pl.getExpectedAverageDegree()))


			gen.generatePowerlawCommunitySizeSequence(minCommunitySize=communityMinSize, maxCommunitySize=communityMaxSize, communitySizeExp=-1 * communityExp)
			# mixing parameter
			#print("mixing parameter")
			gen.setMu(mu)
			# Add some small constants to the parameters for the reference implementation to
			# ensure it won't say the average degree is too low.
			refImplParams = "-N {0} -k {1} -maxk {2} -mu {3} -minc {4} -maxc {5} -t1 {6} -t2 {7}".format(n * scale, avgDegree + 1e-4, maxDegree, mu, max(communityMinSize, 3), communityMaxSize, nodeDegreeExp + 0.001, communityExp)
			cls.params["refImplParams"] = refImplParams
			print(refImplParams)
		else:
			if scale > 1:
				# scale communities
				cData = communities.getVector()
				cDataCopy = cData[:]
				b = communities.upperBound()
				for s in range(1, scale):
					cDataExtend = [i + (b * s) for i in cDataCopy]
					cData = cData + cDataExtend
				assert (len(cData) == n * scale)
				gen.setPartition(Partition(0, cData))
			else:
				gen.setPartition(communities)
			# degree sequence
			gen.setDegreeSequence(degSeq * scale)
			# mixing parameter
			localCoverage = LocalPartitionCoverage(G, communities).run().scores()
			gen.setMu([1.0 - x for x in localCoverage] * scale)
		return gen

cdef extern from "<networkit/generators/MocnikGenerator.hpp>":

	cdef cppclass _MocnikGenerator "NetworKit::MocnikGenerator"(_StaticGraphGenerator):
		_MocnikGenerator(count dim, count n, double k, bool_t weighted) except +
		_MocnikGenerator(count dim, vector[count] ns, double k, bool_t weighted) except +
		_MocnikGenerator(count dim, vector[count] ns, vector[double] ks, bool_t weighted) except +
		_MocnikGenerator(count dim, count n, double k, vector[double] weighted) except +
		_MocnikGenerator(count dim, vector[count] ns, double k, vector[double] weighted) except +
		_MocnikGenerator(count dim, vector[count] ns, vector[double] ks, vector[double] weighted) except +

cdef class MocnikGenerator(StaticGraphGenerator):
	"""
	MocnikGenerator(dim, n, k, weighted)

	Creates random spatial graphs according to the Mocnik model (improved algorithm).

	Please cite the following publications, in which you will find a
	description of the model:

	Franz-Benjamin Mocnik: "The Polynomial Volume Law of Complex Networks in
	the Context of Local and Global Optimization", Scientific Reports 8(11274)
	2018. doi: 10.1038/s41598-018-29131-0

	Franz-Benjamin Mocnik, Andrew Frank: "Modelling Spatial Structures",
	Proceedings of the 12th Conference on Spatial Information Theory (COSIT),
	2015, pages 44-64. doi: 10.1007/978-3-319-23374-1_3

	Parameters
	----------
	dim : int
		Dimension of the space.
	n : int
		Number of nodes in the graph; or a list containing the numbers
		of nodes in each layer in case of a hierarchical model.
	k : float
		Density parameter, determining the ratio of edges to nodes; in
		case of a hierarchical model, also a list of density parameters can be
		provided.
 	weighted : bool, optional
		Determines whether weights should be added to the edges;
		in case of a hierarchical model, also a list of relative weights can be
		provided. Default: False
	"""

	def __cinit__(self, dim, n, k, weighted=False):
		if dim < 1:
			raise ValueError("Dimension must be > 0")
		elif (type(n) is int) and (type(k) is float or type(k) is int) and (weighted is False or weighted is True):
			self._this = new _MocnikGenerator(<count> dim, <count> n, <double> k, <bool_t> weighted)
		elif (type(n) is list) and all(type(item) is int for item in n) and (type(k) is float or type(k) is int) and (weighted is False or weighted is True):
			self._this = new _MocnikGenerator(<count> dim, <vector[count]> n, <double> k, <bool_t> weighted)
		elif (type(n) is list) and all(type(item) is int for item in n) and (type(k) is list) and all(type(item) is float or type(item) is int for item in k) and (weighted is False or weighted is True):
			self._this = new _MocnikGenerator(<count> dim, <vector[count]> n, <vector[double]> k, <bool_t> weighted)
		elif (type(n) is int) and (type(k) is float or type(k) is int) and (type(weighted) is list) and all(type(item) is float or type(item) is int for item in weighted):
			self._this = new _MocnikGenerator(<count> dim, <count> n, <double> k, <vector[double]> weighted)
		elif (type(n) is list) and all(type(item) is int for item in n) and (type(k) is float or type(k) is int) and (type(weighted) is list) and all(type(item) is float or type(item) is int for item in weighted):
			self._this = new _MocnikGenerator(<count> dim, <vector[count]> n, <double> k, <vector[double]> weighted)
		elif (type(n) is list) and all(type(item) is int for item in n) and (type(k) is list) and all(type(item) is float or type(item) is int for item in k) and (type(weighted) is list) and all(type(item) is float or type(item) is int for item in weighted):
			self._this = new _MocnikGenerator(<count> dim, <vector[count]> n, <vector[double]> k, <vector[double]> weighted)
		else:
			pass

cdef extern from "<networkit/generators/MocnikGeneratorBasic.hpp>":

	cdef cppclass _MocnikGeneratorBasic "NetworKit::MocnikGeneratorBasic"(_StaticGraphGenerator):
		_MocnikGeneratorBasic(count dim, count n, double k) except +

cdef class MocnikGeneratorBasic(StaticGraphGenerator):
	"""
	MocnikGeneratorBasic(dim, n, k)

	Creates random spatial graphs according to the Mocnik model (non-improved algorithm).

	Please cite the following publications, in which you will find a
	description of the model:

	Franz-Benjamin Mocnik: "The Polynomial Volume Law of Complex Networks in
	the Context of Local and Global Optimization", Scientific Reports 8(11274)
	2018. doi: 10.1038/s41598-018-29131-0

	Franz-Benjamin Mocnik, Andrew Frank: "Modelling Spatial Structures",
	Proceedings of the 12th Conference on Spatial Information Theory (COSIT),
	2015, pages 44-64. doi: 10.1007/978-3-319-23374-1_3

	Parameters
	----------
	dim : int
		Dimension of the space.
	n : int
		Number of nodes in the graph.
	k : float
		Density parameter, determining the ratio of edges to nodes.
	"""

	def __cinit__(self, dim, n, k):
		self._this = new _MocnikGeneratorBasic(dim, n, k)

cdef extern from "<networkit/generators/HavelHakimiGenerator.hpp>":

	cdef cppclass _HavelHakimiGenerator "NetworKit::HavelHakimiGenerator"(_StaticGraphGenerator):
		_HavelHakimiGenerator(vector[count] degreeSequence, bool_t ignoreIfRealizable) except +
		bool_t isRealizable() except +
		bool_t getRealizable() except +

cdef class HavelHakimiGenerator(StaticGraphGenerator):
	""" 
	HavelHakimiGenerator(sequence, ignoreIfRealizable=True)
	
	Havel-Hakimi algorithm for generating a graph according to a given degree sequence.

	The sequence, if it is realizable, is reconstructed exactly. The resulting graph usually
	has a high clustering coefficient. Construction runs in linear time O(m).

	If the sequence is not realizable, depending on the parameter ignoreIfRealizable, either
	an exception is thrown during generation or the graph is generated with a modified degree
	sequence, i.e. not all nodes might have as many neighbors as requested.

	Parameters
	----------
	sequence : list(int)
		Degree sequence to realize. Must be non-increasing.
	ignoreIfRealizable : bool, optional
		If True, generate the graph even if the degree sequence is not realizable. Some nodes may get lower degrees than requested in the sequence. Default: True
	"""

	def __cinit__(self, vector[count] degreeSequence, ignoreIfRealizable=True):
		self._this = new _HavelHakimiGenerator(degreeSequence, ignoreIfRealizable)

	def isRealizable(self):
		"""
		isRealizable()

		Test if degree sequence is realizable.

		Returns
		-------
		bool
			Indicator for realizable degree sequence.
		"""
		return (<_HavelHakimiGenerator*>(self._this)).isRealizable()

	def getRealizable(self):
		"""
		getRealizable()

		Get realizable state without testing.

		Returns
		-------
		bool
			Indicator for realizable degree sequence.
		"""
		return (<_HavelHakimiGenerator*>(self._this)).getRealizable()

	@classmethod
	def fit(cls, Graph G, scale=1):
		degSeq = DegreeCentrality(G).run().scores()
		return cls(degSeq * scale, ignoreIfRealizable=True)

cdef extern from "<networkit/generators/DynamicHyperbolicGenerator.hpp>":

	cdef cppclass _DynamicHyperbolicGenerator "NetworKit::DynamicHyperbolicGenerator":
		_DynamicHyperbolicGenerator(count numNodes, double avgDegree, double gamma, double T, double moveEachStep, double moveDistance) except +
		vector[_GraphEvent] generate(count nSteps) except +
		_Graph getGraph() except +
		vector[_Point2D] getCoordinates() except +

cdef class DynamicHyperbolicGenerator:
	""" 
	DynamicHyperbolicGenerator(numNodes, avgDegree=6.0, gamma=3.0, T=0.0, moveEachStep=1.0, moveDistance=0.1)

	Dynamic graph generator according to the hyperbolic unit disk model.

	Parameters
	----------
	numNodes : int
		Number of nodes.
	avgDegree : float
		Average degree of the resulting graph. Default: 6.0
	gamma : float
		Power-law exponent of the resulting graph. Default: 3.0
	T : float
		Temperature, selecting a graph family on the continuum between hyperbolic unit disk graphs and Erdos-Renyi graphs. Default: 0.0
	moveEachStep : float
		Fraction of nodes to be moved in each time step. The nodes are chosen randomly each step. Default: 1.0
	moveDistance: float
		Base value for the node movements. Default: 0.1
	"""
	cdef _DynamicHyperbolicGenerator* _this

	def __cinit__(self, numNodes, avgDegree = 6, gamma = 3, T = 0, moveEachStep = 1, moveDistance = 0.1):
		if gamma <= 2:
				raise ValueError("Exponent of power-law degree distribution must be > 2")
		self._this = new _DynamicHyperbolicGenerator(numNodes, avgDegree = 6, gamma = 3, T = 0, moveEachStep = 1, moveDistance = 0.1)

	def __dealloc__(self):
		del self._this

	def generate(self, nSteps):
		""" 
		generate(nSteps)
		
		Generate event stream.

		Parameters
		----------
		nSteps : int
			Number of time steps in the event stream.

		Returns
		-------
		list(networkit.dynamics.GraphEvent)
			List of graph events.
		"""
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]

	def getGraph(self):
		"""
		getGraph()

		Return current graph.

		Returns
		-------
		networkit.Graph
			The current graph.
		"""
		return Graph().setThis(self._this.getGraph())

	def getCoordinates(self):
		""" 
		getCoordinates()
		
		Get coordinates in the Poincare disk.
		
		Returns
		-------
		list(tuple(float,float))
			2D coordinates for every node in the graph.
		"""
		return toPoint2DVector(self._this.getCoordinates())

cdef extern from "<networkit/generators/DynamicDorogovtsevMendesGenerator.hpp>":

	cdef cppclass _DynamicDorogovtsevMendesGenerator "NetworKit::DynamicDorogovtsevMendesGenerator":
		_DynamicDorogovtsevMendesGenerator() except +
		vector[_GraphEvent] generate(count nSteps) except +


cdef class DynamicDorogovtsevMendesGenerator:
	""" 
	DynamicDorogovtsevMendesGenerator()

	Generates a graph according to the Dorogovtsev-Mendes model.
	"""
	cdef _DynamicDorogovtsevMendesGenerator* _this

	def __cinit__(self):
		self._this = new _DynamicDorogovtsevMendesGenerator()

	def __dealloc__(self):
		del self._this

	def generate(self, nSteps):
		""" 
		generate(nSteps)
		
		Generate event stream.

		Parameters
		----------
		nSteps : int
			Number of time steps in the event stream.

		Returns
		-------
		list(networkit.dynamics.GraphEvent)
			List of graph events.
		"""
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]

cdef extern from "<networkit/generators/RmatGenerator.hpp>":

	cdef cppclass _RmatGenerator "NetworKit::RmatGenerator"(_StaticGraphGenerator):
		_RmatGenerator(count scale, count edgeFactor, double a, double b, double c, double d, bool_t weighted, count reduceNodes) except +

cdef class RmatGenerator(StaticGraphGenerator):
	"""
	RmatGenerator(scale, edgeFactor, a, b, c, d, weighted=False, reduceNodes=0)

	Generates static R-MAT graphs. R-MAT (recursive matrix) graphs are
	random graphs with n=2^scale nodes and m=nedgeFactor edges.
	More details at http://www.graph500.org or in the original paper:
	Deepayan Chakrabarti, Yiping Zhan, Christos Faloutsos:
	R-MAT: A Recursive Model for Graph Mining. SDM 2004: 442-446.

	Parameters
	----------
	scale : int
		Number of nodes = 2^scale
	edgeFactor : int
		Number of edges = number of nodes * edgeFactor
	a : float
		Probability for quadrant upper left
	b : float
		Probability for quadrant upper right
	c : float
		Probability for quadrant lower left
	d : float
		Probability for quadrant lower right
	weighted : bool, optional
		Indicates whether the resulting graph should be weighted. Default: False
	reduceNodes : int, optional
		The number of nodes, which should be deleted from the generated graph. Default: 0
	"""
	paths = {"kronfitPath" : None}

	def __cinit__(self, count scale, count edgeFactor, double a, double b, double c, double d, bool_t weighted=False, count reduceNodes=0):
		self._this = new _RmatGenerator(scale, edgeFactor, a, b, c, d, weighted, reduceNodes)

	@classmethod
	def setPaths(cls, kronfitPath):
		"""
		setPaths(kronfitPath)

		Helper function for fit. Sets an output folder for writing the result to.

		Parameters
		----------
		kronfitPath : str
			String containing a path.
		"""
		cls.paths["kronfitPath"] = kronfitPath

	@classmethod
	def fit(cls, G, scale=1, initiator=None, kronfit=True, iterations=50):
		"""
		fit(G, scale=1, initiator=None, kronfit=True, iterations=50)

		Fit model to input graph

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		scale : int, optional
			Scale the maximum number of nodes by a factor. Default: 1
		initiator : tuple(float, float, float, float), optional
			Initiate quadrants with custom values. Default: None
		kronfit : bool, optional
			Indicates whether a slower but more accurate fitting functions is used. Default: True
		iterations: int, optional
			Number of iterations. Default: 50

		Returns
		-------
		networkit.Graph
			A new scaled graph.
		"""
		import math
		import re
		import subprocess
		import os
		import random
		from networkit import graphio
		if initiator:
			(a,b,c,d) = initiator
		else:
			if kronfit:
				with tempfile.TemporaryDirectory() as tmpdir:
					if cls.paths["kronfitPath"] is None:
						raise RuntimeError("call setPaths class method first to configure")
					# write graph
					tmpGraphPath = os.path.join(tmpdir, "tmp.edgelist")
					tmpOutputPath = os.path.join(tmpdir, "tmp.kronfit")
					graphio.writeGraph(G, tmpGraphPath, graphio.Format.EdgeList, separator="\t", firstNode=1, bothDirections=True)
					# call kronfit
					args = [cls.paths["kronfitPath"], "-i:{0}".format(tmpGraphPath), "-gi:{0}".format(str(iterations)), "-o:{}".format(tmpOutputPath)]
					subprocess.call(args)
					# read estimated parameters
					with open(tmpOutputPath) as resultFile:
						for line in resultFile:
							if "initiator" in line:
								matches = re.findall("\d+\.\d+", line)
								weights = [float(s) for s in matches]
			else:
				# random weights because kronfit is slow
				weights = (random.random(), random.random(), random.random(), random.random())
			# normalize
			nweights = [w / sum(weights) for w in weights]
			(a,b,c,d) = nweights
		print("using initiator matrix [{0},{1};{2},{3}]".format(a,b,c,d))
		# other parameters
		(n,m) = GraphTools.size(G)
		scaleParameter = math.ceil(math.log(n * scale, 2))
		edgeFactor = math.floor(m / n)
		reduceNodes = (2**scaleParameter) - (scale * n)
		print("random nodes to delete to achieve target node count: ", reduceNodes)
		return RmatGenerator(scaleParameter, edgeFactor, a, b, c, d, False, reduceNodes)

cdef extern from "<networkit/generators/DynamicForestFireGenerator.hpp>":

	cdef cppclass _DynamicForestFireGenerator "NetworKit::DynamicForestFireGenerator":
		_DynamicForestFireGenerator(double p, bool_t directed, double r) except +
		vector[_GraphEvent] generate(count nSteps) except +
		_Graph getGraph() except +


cdef class DynamicForestFireGenerator:
	""" 
	DynamicForestFireGenerator(p, directed, r = 1.0)

	Generates a graph according to the forest fire model. The forest fire generative model produces dynamic 
	graphs with the properties heavy tailed, degree distribution communities, densification, power law, shrinking diameter.

	See Leskovec, Kleinberg, Faloutsos: Graphs over Tim: Densification Laws,
	Shringking Diameters and Possible Explanations


	Parameters
	----------
	p : float
		Forward burning probability.
	directed : bool
		Decides whether the resulting graph should be directed.
	r : float, optional
		Backward burning probability. Default 1.0
	"""
	cdef _DynamicForestFireGenerator* _this

	def __cinit__(self, p, directed, r = 1.0):
		self._this = new _DynamicForestFireGenerator(p, directed, r)

	def __dealloc__(self):
		del self._this

	def generate(self, nSteps):
		"""
		generate(nSteps)

		Generate event stream.

		Parameters
		----------
		nSteps : int
			Number of time steps in the event stream.
		
		Returns
		-------
		list(networkit.dynamics.GraphEvent)
			List of graph events.
		"""
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]

cdef extern from "<networkit/generators/RegularRingLatticeGenerator.hpp>":

	cdef cppclass _RegularRingLatticeGenerator "NetworKit::RegularRingLatticeGenerator"(_StaticGraphGenerator):
		_RegularRingLatticeGenerator(count nNodes, count nNeighbors) except +

cdef class RegularRingLatticeGenerator(StaticGraphGenerator):
	"""
	RegularRingLatticeGenerator(nNodes, nNeighbors)

	Constructs a regular ring lattice.

	Parameters
	----------
	nNodes : int
		Number of nodes in the target graph.
	nNeighbors : int
		Number of neighbors on each side of a node.
	"""

	def __cinit__(self, nNodes, nNeighbors):
		self._this = new _RegularRingLatticeGenerator(nNodes, nNeighbors)


cdef extern from "<networkit/generators/WattsStrogatzGenerator.hpp>":

	cdef cppclass _WattsStrogatzGenerator "NetworKit::WattsStrogatzGenerator"(_StaticGraphGenerator):
		_WattsStrogatzGenerator(count nNodes, count nNeighbors, double p) except +

cdef class WattsStrogatzGenerator(StaticGraphGenerator):
	""" 
	WattsStrogatzGenerator(nNodes, nNeighbors, p)

	Generates a graph according to the Watts-Strogatz model.

	First, a regular ring lattice is generated. Then edges are rewired
	with a given probability.

	Parameters
	----------
	nNodes : int
		Number of nodes in the target graph.
	nNeighbors : int
		Number of neighbors on each side of a node.
	p : float
		Rewiring probability.
	"""

	def __cinit__(self, nNodes, nNeighbors, p):
		self._this = new _WattsStrogatzGenerator(nNodes, nNeighbors, p)

cdef extern from "<networkit/generators/EdgeSwitchingMarkovChainGenerator.hpp>":

	cdef cppclass _EdgeSwitchingMarkovChainGenerator "NetworKit::EdgeSwitchingMarkovChainGenerator"(_StaticGraphGenerator):
		_EdgeSwitchingMarkovChainGenerator(vector[count] degreeSequence, bool_t ignoreIfNotRealizable, count numSwitches) except +
		bool_t isRealizable() except +
		bool_t getRealizable() except +

cdef class EdgeSwitchingMarkovChainGenerator(StaticGraphGenerator):
	"""
	EdgeSwitchingMarkovChainGenerator(degreeSequence, ignoreIfNotRealizable=False, numSwitchesPerEdge=10)

	Graph generator for generating a random simple graph with exactly the given degree sequence based on the Edge-Switching Markov-Chain method.

	This implementation is based on the paper
	"Random generation of large connected simple graphs with prescribed degree distribution" by Fabien Viger and Matthieu Latapy,
	available at http://www-rp.lip6.fr/~latapy/FV/generation.html, however without preserving connectivity (this could later be added as
	optional feature).

	The Havel-Hakami generator is used for the initial graph generation, then the Markov-Chain Monte-Carlo algorithm as described and
	implemented by Fabien Viger and Matthieu Latapy but without the steps for ensuring connectivity is executed. This should lead to a
	graph that is drawn uniformly at random from all graphs with the given degree sequence.

	Note that at most 10 times the number of edges edge swaps are performed (same number as in the abovementioned implementation) and
	in order to limit the running time, at most 200 times as many attempts to perform an edge swap are made (as certain degree distributions
	do not allow edge swaps at all).

	Parameters
	----------
	degreeSequence : list(int)
		The degree sequence that shall be generated.
	ignoreIfNotRealizable : bool, optional
		If true, generate the graph even if the degree sequence is not realizable. 
		Some nodes may get lower degrees than requested in the sequence. Default: False
	numSwitchesPerEdge : int, optional
		Average number of edge switches per edge produced. Default: 10
	"""

	def __cinit__(self, vector[count] degreeSequence, bool_t ignoreIfNotRealizable = False, count numSwitchesPerEdge = 10):
		self._this = new _EdgeSwitchingMarkovChainGenerator(degreeSequence, ignoreIfNotRealizable, numSwitchesPerEdge)

	def isRealizable(self):
		"""
		isRealizable()

		Test if degree sequence is realizable.

		Returns
		-------
		bool
			Indicator for realizable degree sequence.
		"""
		return (<_EdgeSwitchingMarkovChainGenerator*>(self._this)).isRealizable()

	def getRealizable(self):
		"""
		getRealizable()

		Get realizable state without testing.

		Returns
		-------
		bool
			Indicator for realizable degree sequence.
		"""
		return (<_EdgeSwitchingMarkovChainGenerator*>(self._this)).getRealizable()

	@classmethod
	def fit(cls, Graph G, scale=1):
		"""
		fit(G, scale=1)

		Fit model to input graph

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		scale : int, optional
			Scale the maximum number of nodes by a factor. Default: 1

		Returns
		-------
		networkit.Graph
			A new scaled graph.
		"""
		degSeq = DegreeCentrality(G).run().scores()
		return cls(degSeq * scale, ignoreIfRealizable=True)

ConfigurationModelGenerator = EdgeSwitchingMarkovChainGenerator

cdef extern from "<networkit/generators/DynamicPathGenerator.hpp>":

	cdef cppclass _DynamicPathGenerator "NetworKit::DynamicPathGenerator":
		_DynamicPathGenerator() except +
		vector[_GraphEvent] generate(count nSteps) except +

cdef class DynamicPathGenerator:
	""" 
	DynamicPathGenerator()

	Example dynamic graph generator: Generates a dynamically growing path. 
	"""
	cdef _DynamicPathGenerator* _this

	def __cinit__(self):
		self._this = new _DynamicPathGenerator()

	def __dealloc__(self):
		del self._this

	def generate(self, nSteps):
		"""
		generate(nSteps)

		Generate event stream.

		Parameters
		----------
		nSteps : int
			Number of time steps in the event stream.
		
		Returns
		-------
		list(networkit.dynamics.GraphEvent)
			List of graph events.
		"""
		return [GraphEvent(ev.type, ev.u, ev.v, ev.w) for ev in self._this.generate(nSteps)]

class BTERReplicator:
	"""
	BTERReplicator()

	Wrapper class that calls the BTER graph generator implementation in
	FEASTPACK from http://www.sandia.gov/~tgkolda/feastpack/ using GNU
	Octave.

	Note that BTER needs the rng method which is unavailable in Octave, but
	the call in bter.m can be easily replaced.
	"""
	matlabname = 'octave'
	matlabScript = """
	addpath('{0}');
	filename = 'bter_input.mat';
	load(filename);
	addpath('{1}');
	tStart = tic;
	[ccd,gcc] = ccperdeg(G);
	nd = accumarray(nonzeros(sum(G,2)),1);
	nd = nd * {2};
	tFit = toc(tStart);
	tStart = tic;
	[E1,E2] = bter(nd,ccd,'verbose',false,'blowup',10);
	tGenerate = toc(tStart);
	G_bter = bter_edges2graph(E1,E2);
	save('-v7', '{3}', 'G_bter', 'tFit', 'tGenerate');
	exit;
	"""
	feastpackPath = "."


	@classmethod
	def setPaths(cls, feastpackPath):
		"""
		setPaths(feastpackPath)

		Helper function for fit. Sets an output folder for writing the result to.

		Parameters
		----------
		feastpackPath : str
			String containing a path.
		"""
		cls.feastpackPath = feastpackPath

	def __init__(self, G, scale=1):
		self.G = G
		self.scale = scale

	def generate(self):
		"""
		generate()

		Generate graph.
		"""
		from . import graphio
		with tempfile.TemporaryDirectory() as tmpdir:
			scriptPath = os.path.join(tmpdir, "bter_wrapper.m")
			tempFileOut = os.path.join(tmpdir, 'bter_output.mat')
			tempFileIn = os.path.join(tmpdir, 'bter_input.mat')
			# write MATLAB script
			with open(scriptPath, 'w') as matlabScriptFile:
				matlabScriptFile.write(self.matlabScript.format(tmpdir, self.feastpackPath, self.scale, tempFileOut))
			graphio.writeMat(self.G, tempFileIn)
			subprocess.call([self.matlabname, '-qf', scriptPath])
			G_bter = graphio.readMat(tempFileOut, key='G_bter')
			matlabObject = scipy.io.loadmat(tempFileOut)
			self.t_fit = matlabObject["tFit"][0][0]
			self.t_generate = matlabObject["tGenerate"][0][0]
			return G_bter

	@classmethod
	def fit(cls, G, scale=1):
		"""
		fit(G, scale=1)

		Fit model to input graph

		Parameters
		----------
		G : networkit.Graph
			The input graph.
		scale : int, optional
			Scale the maximum number of nodes by a factor. Default: 1

		Returns
		-------
		networkit.Graph
			A new scaled graph.
		"""
		return cls(G, scale)
