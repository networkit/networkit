# distutils: language=c++

from libc.stdint cimport uint64_t
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph

cdef extern from "<networkit/randomization/GlobalCurveball.hpp>":

	cdef cppclass _GlobalCurveball "NetworKit::GlobalCurveball"(_Algorithm):
		_GlobalCurveball(_Graph, count, bool_t, bool_t) except +
		_Graph getGraph() except +

cdef class GlobalCurveball(Algorithm):
	"""
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

	Parameters:
	-----------

	G : networkit.Graph
		The graph to be randomized. For a given degree sequence, e.g.
		generators.HavelHakimi can be used to obtain this graph.

	number_of_global_rounds:
		Number of global rounds to carry out. The runtime scales
		asymptotically linearly in this parameter. Default: 20,
		which yields good results experimentally (see Paper).

	allowSelfLoops:
		Has to be False for undirected graphs. For directed graphs
		the randomization Markov chain is only irreducible if self loops
		are allows. If they are forbidden, the degreePreservingShuffle
		proprocessing has to be enabled. Otherwhise, not all topologies
		can be produced.

	degreePreservingShufflePreprocessing:
		Execute the DegreePreservingShuffle algorithm before executing
		Global Curveball. It's more efficient than manually invoking
		the algorithm.

	Warning:
	--------
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

	"""

	Get randomized graph after invocation of run().

	"""

	def getGraph(self):
		return Graph().setThis((<_GlobalCurveball*>self._this).getGraph())

cdef extern from "<networkit/randomization/CurveballUniformTradeGenerator.hpp>":

	cdef cppclass _CurveballUniformTradeGenerator "NetworKit::CurveballUniformTradeGenerator":
		_CurveballUniformTradeGenerator(count runLength, count numNodes) except +
		vector[pair[node, node]] generate() nogil except +

cdef class CurveballUniformTradeGenerator:

	"""
	Generates a trade sequence consisting of num_trades many single trades.
	Each trade contains two different node indices drawn uniformly at random
	from the interval [0, num_nodes).

	Parameters:
	-----------

	num_trades:
	   Number of trades to generate.

	num_nodes:
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
	Generates a trade sequence consisting of num_global_trades global trades
	targeting node ids from the range [0, num_nods).

	If you are only using this generator, consider using the GlobalCurveball
	algorithm directly as it has a better performance / memory footprint.

	Parameters:
	-----------

	num_global_trades:
	   Number of global trades to generate (i.e. the resulting sequence contains
	   num_global_trades * floor(num_nodes / 2) trades)

	num_nodes:
	   Number of node indices to draw from

	"""
	cdef _CurveballGlobalTradeGenerator *_this

	def __cinit__(self, count num_global_trades, count num_nodes):
		self._this = new _CurveballGlobalTradeGenerator(num_global_trades, num_nodes)

	def __dealloc__(self):
		del self._this

	def generate(self):
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

	Parameters:
	-----------

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
		with nogil:
			(<_Curveball*>(self._this)).run(trades)
		return self

	def getGraph(self):
		return Graph().setThis((<_Curveball*>self._this).getGraph())

	def getNumberOfAffectedEdges(self):
		return (<_Curveball*>(self._this)).getNumberOfAffectedEdges()

cdef extern from "<networkit/randomization/DegreePreservingShuffle.hpp>":
	cdef cppclass _DegreePreservingShuffle "NetworKit::DegreePreservingShuffle"(_Algorithm):
		_DegreePreservingShuffle(_Graph) except +
		_Graph getGraph() except +
		vector[node] getPermutation() except +

cdef class DegreePreservingShuffle(Algorithm):
	"""
	Implementation of the preprocessing step proposed in
	"Smaller Universes for Uniform Sampling of 0,1-matrices with fixed row and column sums"
	by Annabell Berger, Corrie Jacobien Carstens [https://arxiv.org/abs/1803.02624]

	The algorithms randomizes a graph without changing its topology simply
	by renaming nodes. For any degree d (in case of an directed graph it's a degree pair)
	consider the set X_d of node ids which have this degree. Then shuffle the ids in X_d.

	Hence the algorithm satisfies: For all x in Ginput:
	 i)  Ginput.degreeIn(x) = Goutput.degreeIn(x)
	 ii) Ginput.degreeOut(x) = Goutput.degreeOut(x)

	The authors argue that applying this preprocessing step before executing (Global)Curveball
	leads to a faster mixing time. If you want to use it as a preprocessing step to GlobalCurveball,
	it's more efficient to set degreePreservingShufflePreprocessing in GlobalCurveball's constructor.

	Parameters:
	-----------

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
		return Graph().setThis((<_DegreePreservingShuffle*>self._this).getGraph())

	def getPermutation(self):
		return (<_DegreePreservingShuffle*>(self._this)).getPermutation()
