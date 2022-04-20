# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport count, index, node

cdef extern from "<networkit/randomization/EdgeSwitching.hpp>":
	cdef cppclass _EdgeSwitching "NetworKit::EdgeSwitching"(_Algorithm):
		_EdgeSwitching(_Graph, double, bool_t) except +
		void run(count) nogil except +
		_Graph getGraph() except +
		count getNumberOfAffectedEdges()
		double getNumberOfSwitchesPerEdge()
		void setNumberOfSwitchesPerEdge(double)

	cdef cppclass _EdgeSwitchingInPlace "NetworKit::EdgeSwitchingInPlace"(_Algorithm):
		_EdgeSwitchingInPlace(_Graph, double) except +
		void run(count) nogil except +
		count getNumberOfAffectedEdges()
		double getNumberOfSwitchesPerEdge()
		void setNumberOfSwitchesPerEdge(double)


cdef class EdgeSwitching(Algorithm):
	"""
	EdgeSwitching(G, numberOfSwapsPerEdge=10.0, degreePreservingShufflePreprocessing=True)

	The Edge Switching Markov Chain ["The markov chain simulation method for generating connected
	power law random graphs", Mihail and Zegura] perturbs simple directed or undirected graphs
	while preserving their degrees. In each step, we select two edges uniformly at random, and
	exchange their endpoints. Swaps that introduce multi-edges or self-loops are rejected WITHOUT
	replacement -- this is necessary to allow uniform sampling [see "Switching edges to randomize
	networks: what goes wrong and how to fix it", Carstens and Horadam]. The number of successful
	swaps can be queried using getNumberOfAffectedEdges()/2.

	We provide two implementations: EdgeSwitching takes a copy of the input graph and is more
	versatile; EdgeSwitchingInPlace works directly on the graph supplied but cannot carry-out
	an initial degreePreservingShufflePreprocessing.

	In general, simple edge switching does not yield a uniform distribution on simple DIRECTED
	graphs because the orientation of directed triangles cannot be changed. Using
	DegreePreservingShuffle as a preprocessing step overcomes this limitation. The
	preprocessing can also jump-start the perturbation process, yielding to potentially faster
	mixing. It is performed by default for owned graphs.

	Parameters
	----------
	G : networkit.Graph
		The graph to be randomized.
	numberOfSwapsPerEdge : float, optional
		The average number of swaps to be carried out per edge.
		Has to be non-negative. Default: 10.0
	degreePreservingShufflePreprocessing : bool, optional
		If true (default), in a preprocessing step jump starts the perturbation process.
		For undirected graph, this yields faster mixing; for directed graphs, it is
		necessary in order to obtain an unbiased sampling. Default: True
	"""

	def __cinit__(self, G, numberOfSwapsPerEdge = 10.0, degreePreservingShufflePreprocessing = True):
		if isinstance(G, Graph):
			self._this = new _EdgeSwitching((<Graph>G)._this, numberOfSwapsPerEdge, degreePreservingShufflePreprocessing)
		else:
			raise RuntimeError("Parameter G has to be a graph")

	def run(self):
		"""Perform edge switching. May be called multiple times."""
		(<_EdgeSwitching*>self._this).run()

	def getGraph(self):
		return Graph().setThis((<_EdgeSwitching*>self._this).getGraph())

	def getNumberOfAffectedEdges(self):
		return (<_EdgeSwitching*>(self._this)).getNumberOfAffectedEdges()

	def getNumberOfSwitchesPerEdge(self):
		return (<_EdgeSwitching*>(self._this)).getNumberOfSwitchesPerEdge()

	def setNumberOfSwitchesPerEdge(self, numberOfSwitchesPerEdge):
		(<_EdgeSwitching*>(self._this)).setNumberOfSwitchesPerEdge(numberOfSwitchesPerEdge)

cdef class EdgeSwitchingInPlace(Algorithm):
	"""
	EdgeSwitchingInPlace(G, numberOfSwitchesPerEdge)
	The Edge Switching Markov Chain ["The markov chain simulation method for generating connected
	power law random graphs", Mihail and Zegura] perturbs simple directed or undirected graphs
	while preserving their degrees. In each step, we select two edges uniformly at random, and
	exchange their endpoints. Swaps that introduce multi-edges or self-loops are rejected WITHOUT
	replacement -- this is necessary to allow uniform sampling [see "Switching edges to randomize
	networks: what goes wrong and how to fix it", Carstens and Horadam]. The number of successful
	swaps can be queried using getNumberOfAffectedEdges()/2.

	We provide two implementations: EdgeSwitching takes a copy of the input graph and is more
	versatile; EdgeSwitchingInPlace works directly on the graph supplied but cannot carry-out
	an initial degreePreservingShufflePreprocessing.

	In general, simple edge switching does not yield a uniform distribution on simple DIRECTED
	graphs because the orientation of directed triangles cannot be changed. Using
	DegreePreservingShuffle as a preprocessing step overcomes this limitation. The
	preprocessing can also jump-start the perturbation process, yielding to potentially faster
	mixing. It is only available for EdgeSwitching.

	The implementation keeps a local reference (accessible via getGraph) to prevent premature
	garbage collection.

	Parameters
	----------
	G : networkit.Graph
		The graph to be randomized.
	numberOfSwitchesPerEdge : int, optional
		The average number of switches to be carried out per edge.
		Has to be non-negative. Default: 10
	"""
	cdef Graph _localReference # keep reference counter up to prevent GC of graph

	def __cinit__(self, G, numberOfSwitchesPerEdge = 10.0):
		if isinstance(G, Graph):
			self._this = new _EdgeSwitchingInPlace((<Graph>G)._this, numberOfSwitchesPerEdge)
			self._localReference = G
		else:
			raise RuntimeError("Parameter G has to be a graph")

	def getGraph(self):
		"""
		getGraph()

		Return modified graph.

		Returns
		-------
		networkit.Graph
			The graph after applying edge switches.
		"""
		return self._localReference

	def getNumberOfAffectedEdges(self):
		"""
		getNumberOfAffectedEdges()

		Return number of affected edges.

		Returns
		-------
		int
			Number of edges affected by edge switches.
		"""
		return (<_EdgeSwitchingInPlace*>(self._this)).getNumberOfAffectedEdges()

	def getNumberOfSwitchesPerEdge(self):
		"""
		getNumberOfSwitchesPerEdge()

		Return number of switches per edge.

		Returns
		-------
		list(int)
			Number of switches per edges.
		"""
		return (<_EdgeSwitchingInPlace*>(self._this)).getNumberOfSwitchesPerEdge()

	def setNumberOfSwitchesPerEdge(self, numberOfSwitchesPerEdge):
		(<_EdgeSwitchingInPlace*>(self._this)).setNumberOfSwitchesPerEdge(numberOfSwitchesPerEdge)


cdef extern from "<networkit/randomization/GlobalCurveball.hpp>":

	cdef cppclass _GlobalCurveball "NetworKit::GlobalCurveball"(_Algorithm):
		_GlobalCurveball(_Graph, count, bool_t, bool_t) except +
		_Graph getGraph() except +

cdef class GlobalCurveball(Algorithm):
	"""
	GlobalCurveball(G, number_of_global_rounds=20, allowSelfLoops=False, degreePreservingShufflePreprocessing=True)
	Implementation of EM-GCB proposed in "Parallel and I/O-efficient
	Randomisation of Massive Networks using Global Curveball Trades",
	Carstens et al., ESA 2018.

	The algorithm perturbs an unweighted input graph, by iteratively
	randomizing the neighbourhoods of node pairs. For a large number
	of global trades this process is shown to produce an uniform sample
	from the set of all graphs with the same degree sequence as the input
	graph.

	If you do not want to explicitly control the trade sequence,
	we recommend using GlobalCurveball rather than Curveball since
	GlobalCurveball is typically faster and exhibits a smaller memory
	footprint.

	Parameters
	----------

	G : networkit.Graph
		The graph to be randomized. For a given degree sequence, e.g.
		generators.HavelHakimi can be used to obtain this graph.

	number_of_global_rounds : int, optional
		Number of global rounds to carry out. The runtime scales
		asymptotically linearly in this parameter. Default: 20,
		which yields good results experimentally (see Paper).

	allowSelfLoops : bool, optional
		Has to be False for undirected graphs. For directed graphs
		the randomization Markov chain is only irreducible if self loops
		are allows. If they are forbidden, the degreePreservingShuffle
		preprocessing has to be enabled. Otherwise, not all topologies
		can be produced. Default: False

	degreePreservingShufflePreprocessing : bool, optional
		Execute the DegreePreservingShuffle algorithm before executing
		Global Curveball. It's more efficient than manually invoking
		the algorithm. Default: True

	Note
	----
	For directed graphs at least one of allowSelfLoops or
	degreePreservingShufflePreprocessing should be set; for more details
	refer to "Switching edges to randomize networks: what goes wrong
	and how to fix it" by C. J. Carstens K. J. Horadam
	"""
	def __cinit__(self, G, number_of_global_rounds = 20, allowSelfLoops = False, degreePreservingShufflePreprocessing = True):
		if isinstance(G, Graph):
			self._this = new _GlobalCurveball((<Graph>G)._this, number_of_global_rounds, allowSelfLoops, degreePreservingShufflePreprocessing)
		else:
			raise RuntimeError("Parameter G has to be a graph")

	def getGraph(self):
		"""
		getGraph()

		Get randomized graph after invocation of run().

		Returns
		-------
		networkit.Graph
			The randomized graph.
		"""
		return Graph().setThis((<_GlobalCurveball*>self._this).getGraph())

cdef extern from "<networkit/randomization/CurveballUniformTradeGenerator.hpp>":

	cdef cppclass _CurveballUniformTradeGenerator "NetworKit::CurveballUniformTradeGenerator":
		_CurveballUniformTradeGenerator(count runLength, count numNodes) except +
		vector[pair[node, node]] generate() nogil except +

cdef class CurveballUniformTradeGenerator:
	"""
	CurveballUniformTradeGenerator(num_trades, num_nodes)

	Generates a trade sequence consisting of num_trades many single trades.
	Each trade contains two different node indices drawn uniformly at random
	from the interval [0, num_nodes).

	Parameters
	----------
	num_trades : int
	   Number of trades to generate.
	num_nodes : int
	   Number of node indices to draw from
	"""
	cdef _CurveballUniformTradeGenerator *_this

	def __cinit__(self, count num_trades, count num_nodes):
		self._this = new _CurveballUniformTradeGenerator(num_trades, num_nodes)

	def __dealloc__(self):
		del self._this

	def generate(self):
		return self._this.generate()

cdef extern from "<networkit/randomization/CurveballGlobalTradeGenerator.hpp>":

	cdef cppclass _CurveballGlobalTradeGenerator "NetworKit::CurveballGlobalTradeGenerator":
		_CurveballGlobalTradeGenerator(count runLength, count numNodes) except +
		vector[pair[node, node]] generate() nogil except +

cdef class CurveballGlobalTradeGenerator:
	"""
	CurveballGlobalTradeGenerator(num_global_trades, num_nodes)

	Generates a trade sequence consisting of num_global_trades global trades
	targeting node ids from the range [0, num_nods).

	If you are only using this generator, consider using the GlobalCurveball
	algorithm directly as it has a better performance / memory footprint.

	Parameters
	----------
	num_global_trades : int
	   Number of global trades to generate (i.e. the resulting sequence contains
	   num_global_trades * floor(num_nodes / 2) trades)
	num_nodes : int
	   Number of node indices to draw from
	"""
	cdef _CurveballGlobalTradeGenerator *_this

	def __cinit__(self, count num_global_trades, count num_nodes):
		self._this = new _CurveballGlobalTradeGenerator(num_global_trades, num_nodes)

	def __dealloc__(self):
		del self._this

	def generate(self):
		"""
		generate()

		Generate randomized graph.

		Returns
		-------
		networkit.Graph
			The randomized graph.
		"""
		return self._this.generate()

cdef extern from "<networkit/randomization/Curveball.hpp>":

	cdef cppclass _Curveball "NetworKit::Curveball"(_Algorithm):
		_Curveball(_Graph) except +
		void run(vector[pair[node, node]] trades) nogil except +
		_Graph getGraph() except +
		vector[pair[node, node]] getEdges() except +
		count getNumberOfAffectedEdges() except +

cdef class Curveball(Algorithm):
	"""
	Curveball(G)

	Implementation of IM-CB proposed in "Parallel and I/O-efficient
	Randomisation of Massive Networks using Global Curveball Trades",
	Carstens et al., ESA 2018.

	The algorithm perturbs an undirected and unweighted input graph,
	by iteratively randomizing the neighbourhoods of node pairs. For
	a large number of trades this process is shown to produce an
	uniform sample from the set of all graphs with the same degree
	sequence as the input graph.

	If you do not want to explicitly control the trade sequence,
	we recommend using GlobalCurveball rather than Curveball since
	GlobalCurveball is typically faster and exhibits a smaller memory
	footprint.

	Observe that this algorithm does not support the run() method,
	since it requires the trade sequence to be passed. It is possible
	to invoke run(trades) several times, e.g. to reduce the memory
	footprint which increases linearly with the number of trades
	performed in a run.

	Parameters
	----------
	G : networkit.Graph
		The graph to be randomized. For a given degree sequence, e.g.
		generators.HavelHakimi can be used to obtain this graph.
	"""
	def __cinit__(self, G):
		if isinstance(G, Graph):
			self._this = new _Curveball((<Graph>G)._this)
		else:
			raise RuntimeError("Parameter G has to be a graph")

	def run(self, vector[pair[node, node]] trades):
		"""
		run(trades)

		Compute the randomization of the input by given node pairs.

		Parameters
		----------
		trades : list(tuple(int, int))
			List of pairs of nodes used for randomization of the graph.
		"""
		with nogil:
			(<_Curveball*>(self._this)).run(trades)
		return self

	def getGraph(self):
		"""
		getGraph()

		Get randomized graph after invocation of run().

		Returns
		-------
		networkit.Graph
			The randomized graph.
		"""
		return Graph().setThis((<_Curveball*>self._this).getGraph())

	def getNumberOfAffectedEdges(self):
		"""
		getNumberOfAffectedEdges()

		Return number of affected edges.

		Returns
		-------
		int
			Number of edges affected by randomization.
		"""
		return (<_Curveball*>(self._this)).getNumberOfAffectedEdges()

cdef extern from "<networkit/randomization/DegreePreservingShuffle.hpp>":
	cdef cppclass _DegreePreservingShuffle "NetworKit::DegreePreservingShuffle"(_Algorithm):
		_DegreePreservingShuffle(_Graph) except +
		_Graph getGraph() except +
		vector[node] getPermutation() except +

cdef class DegreePreservingShuffle(Algorithm):
	"""
	DegreePreservingShuffle(G)

	Implementation of the preprocessing step proposed in
	"Smaller Universes for Uniform Sampling of 0,1-matrices with fixed row and column sums"
	by Annabell Berger, Corrie Jacobien Carstens [https://arxiv.org/abs/1803.02624]

	The algorithms randomizes a graph without changing its topology simply
	by renaming nodes. For any degree d (in case of an directed graph it's a degree pair)
	consider the set X_d of node ids which have this degree. Then shuffle the ids in X_d.

	Hence the algorithm satisfies: For all x in Ginput:

	1) Ginput.degreeIn(x) = Goutput.degreeIn(x)
	2) Ginput.degreeOut(x) = Goutput.degreeOut(x)

	The authors argue that applying this preprocessing step before executing (Global)Curveball
	leads to a faster mixing time. If you want to use it as a preprocessing step to GlobalCurveball,
	it's more efficient to set degreePreservingShufflePreprocessing in GlobalCurveball's constructor.

	Parameters
	----------
	G : networkit.Graph
		The graph to be randomized. For a given degree sequence, e.g.
		generators.HavelHakimi can be used to obtain this graph.

	"""
	def __cinit__(self, G):
		if isinstance(G, Graph):
			self._this = new _DegreePreservingShuffle((<Graph>G)._this)
		else:
			raise RuntimeError("Parameter G has to be a graph")

	def getGraph(self):
		"""
		getGraph()

		Get randomized graph after invocation of run().

		Returns
		-------
		networkit.Graph
			The randomized graph.
		"""
		return Graph().setThis((<_DegreePreservingShuffle*>self._this).getGraph())

	def getPermutation(self):
		"""
		getPermutation()

		Returns the permutation used for shuffling.

		Returns
		-------
		list(int)
			List of nodes.
		"""
		return (<_DegreePreservingShuffle*>(self._this)).getPermutation()
