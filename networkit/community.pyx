#from .helpers import stdstring

def stdstring(pystring):
	""" convert a Python string to a bytes object which is automatically coerced to std::string"""
	pybytes = pystring.encode("utf-8")
	return pybytes

cdef class CommunityDetector(Algorithm):
	""" Abstract base class for static community detection algorithms """

	def __init__(self, *args, **namedargs):
		if type(self) == CommunityDetector:
			raise RuntimeError("Error, you may not use CommunityDetector directly, use a sub-class instead")

	def getPartition(self):
		"""  Returns a partition of the clustering.

		Returns
		-------
		networkit.Partition:
			A Partition of the clustering.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return Partition().setThis((<_CommunityDetectionAlgorithm*>(self._this)).getPartition())

cdef extern from "<networkit/community/PLM.hpp>":

	cdef cppclass _PLM "NetworKit::PLM"(_CommunityDetectionAlgorithm):
		_PLM(_Graph _G) except +
		_PLM(_Graph _G, bool_t refine, double gamma, string par, count maxIter, bool_t turbo, bool_t recurse) except +
		map[string, vector[count]] getTiming() except +

cdef extern from "<networkit/community/PLM.hpp>" namespace "NetworKit::PLM":

	pair[_Graph, vector[node]] PLM_coarsen "NetworKit::PLM::coarsen" (const _Graph& G, const _Partition& zeta) except +
	_Partition PLM_prolong "NetworKit::PLM::prolong"(const _Graph& Gcoarse, const _Partition& zetaCoarse, const _Graph& Gfine, vector[node] nodeToMetaNode) except +


cdef class PLM(CommunityDetector):
	""" Parallel Louvain Method - the Louvain method, optionally extended to
		a full multi-level algorithm with refinement

		Parameters
		----------
		G : networkit.Graph
			A graph.
		refine : bool, optional
			Add a second move phase to refine the communities.
		gamma : double
			Multi-resolution modularity parameter:
			1.0 -> standard modularity
	 		0.0 -> one community
	 		2m 	-> singleton communities
		par : string
			parallelization strategy
		maxIter : count
			maximum number of iterations for move phase
		turbo : bool, optional
			faster but uses O(n) additional memory per thread
		recurse: bool, optional
			use recursive coarsening, see http://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.049902 for some explanations (default: true)
	"""

	def __cinit__(self, Graph G not None, refine=False, gamma=1.0, par="balanced", maxIter=32, turbo=True, recurse=True):
		self._G = G
		self._this = new _PLM(G._this, refine, gamma, stdstring(par), maxIter, turbo, recurse)

	def getTiming(self):
		"""  Get detailed time measurements.
		"""
		return (<_PLM*>(self._this)).getTiming()

	@staticmethod
	def coarsen(Graph G, Partition zeta, bool_t parallel = False):
		cdef pair[_Graph, vector[node]] result = move(PLM_coarsen(G._this, zeta._this))
		return (Graph().setThis(result.first), result.second)

	@staticmethod
	def prolong(Graph Gcoarse, Partition zetaCoarse, Graph Gfine, vector[node] nodeToMetaNode):
		return Partition().setThis(PLM_prolong(Gcoarse._this, zetaCoarse._this, Gfine._this, nodeToMetaNode))

