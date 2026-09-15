# distutils: language=c++

from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport uint64_t
from cpython.ref cimport PyObject
from cython.operator cimport dereference

from .helpers import stdstring

from .base cimport _Algorithm, Algorithm
from .graph cimport _Graph, Graph
from .structures cimport count, index, node

cdef extern from "cython_helper.h":
	void throw_runtime_error(string message) nogil

cdef cppclass MatchCallbackWrapper:
	void* callback

	__init__(object callback):
		this.callback = <void*>callback

	# This is called within the run() method which is nogil!
	void cython_call_operator(const vector[node]& Match) nogil:
		cdef bool_t error = False
		cdef string message

		# Acquire GIL to allow Python code.
		with gil:
			try:
				(<object>callback)(Match)
			except Exception as e:
				error = True
				message = stdstring(
					"An Exception occurred, aborting execution of iterator: {0}".format(e)
				)

		# Throw only after leaving the with-gil block.
		if error:
			throw_runtime_error(message)


cdef cppclass ParallelMatchCallbackWrapper:
	void* callback

	__init__(object callback):
		this.callback = <void*>callback

	# This is called within the run() method which is nogil!
	void cython_call_operator(index tid, const vector[node]& Match) nogil:
		cdef bool_t error = False
		cdef string message

		# Acquire GIL to allow Python code.
		with gil:
			try:
				(<object>callback)(tid, Match)
			except Exception as e:
				error = True
				message = stdstring(
					"An Exception occurred, aborting execution of iterator: {0}".format(e)
				)

		# Throw only after leaving the with-gil block.
		if error:
			throw_runtime_error(message)

cdef extern from "<networkit/isomorphism/SubgraphIsomorphism.hpp>" namespace "NetworKit":
	cdef enum class _Semantics "NetworKit::SubgraphIsomorphism::Semantics":
		INDUCED,
		MONOMORPHISM

	cdef cppclass _SubgraphIsomorphism "NetworKit::SubgraphIsomorphism" (_Algorithm):
		_SubgraphIsomorphism(const _Graph &pattern, const _Graph &target, _Semantics semantics, count maxMatches) except +
		count numberOfWorkers() except +
		void setNodeLabels(const vector[index] &patternNodeLabels, const vector[index] &targetNodeLabels) except +
		void setEdgeLabels(const vector[index] &patternEdgeLabels, const vector[index] &targetEdgeLabels) except +
		void setCallback(MatchCallbackWrapper callback) except +
		void setCallback(ParallelMatchCallbackWrapper callback) except +
		void setStoreMatches(bool_t storeMatches) except +
		const vector[vector[node]] &getMatches() except +
		count numberOfMatches() except +
		bool_t hasMatch() except +

class Semantics(object):
	"""
	What counts as a match.

	- INDUCED: Pattern edges map to target edges and pattern non-edges map to target non-edges.
	- MONOMORPHISM: Pattern edges map to target edges. Additional target edges are allowed.
	"""
	INDUCED = _Semantics.INDUCED
	MONOMORPHISM = _Semantics.MONOMORPHISM

cdef class SubgraphIsomorphism(Algorithm):
	"""
	Abstract base class for subgraph isomorphism algorithms.

	The algorithms take a small pattern graph and a large target graph and find every match: an
	injective mapping of the pattern nodes to target nodes under which every pattern edge is mapped
	to a target edge. A match is a list indexed by pattern node, so match[u] is the target node that
	pattern node u is mapped to, and it holds networkit.none at ids that are not nodes. Under
	networkit.isomorphism.Semantics.INDUCED, pattern non-edges must also be mapped to target
	non-edges.

	Edge weights are ignored, the pattern must not contain self-loops, and parallel edges are
	collapsed.
	"""

	cdef Graph _pattern
	cdef Graph _target
	cdef MatchCallbackWrapper* _callback
	cdef ParallelMatchCallbackWrapper* _parallelCallback
	cdef object _py_callback

	def __init__(self, *args, **namedargs):
		if type(self) == SubgraphIsomorphism:
			raise RuntimeError("Instantiation of abstract base class")

	def numberOfWorkers(self):
		"""
		numberOfWorkers()

		Returns the number of workers the algorithm uses. The worker id passed to a parallel
		callback is smaller than this number.

		Returns
		-------
		int
			The number of workers; 1 for the sequential algorithms.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).numberOfWorkers()

	def setNodeLabels(self, vector[index] patternNodeLabels, vector[index] targetNodeLabels):
		"""
		setNodeLabels(patternNodeLabels, targetNodeLabels)

		Restricts matches to map every pattern node to a target node with the same label. The label
		networkit.none is a wildcard that matches any label. Passing two empty lists removes the
		node labels. Call this before run().

		Parameters
		----------
		patternNodeLabels : list(int)
			Labels of the pattern nodes, indexed by node id.
		targetNodeLabels : list(int)
			Labels of the target nodes, indexed by node id.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setNodeLabels(patternNodeLabels, targetNodeLabels)

	def setEdgeLabels(self, vector[index] patternEdgeLabels,
					  vector[index] targetEdgeLabels):
		"""
		setEdgeLabels(patternEdgeLabels, targetEdgeLabels)

		Restricts matches to map every pattern edge to a target edge with the same label. The label
		networkit.none is a wildcard that matches any label. Passing two empty lists removes the
		edge labels. Both graphs need edge ids, see networkit.Graph.indexEdges(). Parallel edges
		with different labels are not supported. Call this before run().

		Parameters
		----------
		patternEdgeLabels : list(int)
			Labels of the pattern edges, indexed by edge id.
		targetEdgeLabels : list(int)
			Labels of the target edges, indexed by edge id.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setEdgeLabels(patternEdgeLabels, targetEdgeLabels)

	def setCallback(self, object callback, bool_t parallel=False):
		"""
		setCallback(callback, parallel=False)

		Passes every match to a callback instead of storing it, so getMatches() raises afterwards.
		Replaces a previously set callback. If the callback raises, the search stops and run()
		raises a RuntimeError that carries the original message.

		Parameters
		----------
		callback : callable
			With parallel=False, it is called as callback(match). With parallel=True, it is called
			as callback(workerId, match), where 0 <= workerId < numberOfWorkers().
		parallel : bool, optional
			Whether the callback also receives the worker id. Default: False
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")

		if not callable(callback):
			raise TypeError("callback must be callable")

		# Remove a previously registered wrapper.
		if self._callback != NULL:
			del self._callback
			self._callback = NULL

		if self._parallelCallback != NULL:
			del self._parallelCallback
			self._parallelCallback = NULL

		# Keep the Python callback alive.
		self._py_callback = callback

		if parallel:
			self._parallelCallback = new ParallelMatchCallbackWrapper(callback)
			try:
				(<_SubgraphIsomorphism*>(self._this)).setCallback(
					dereference(self._parallelCallback)
				)
			except BaseException:
				del self._parallelCallback
				self._parallelCallback = NULL
				self._py_callback = None
				raise
		else:
			self._callback = new MatchCallbackWrapper(callback)
			try:
				(<_SubgraphIsomorphism*>(self._this)).setCallback(
					dereference(self._callback)
				)
			except BaseException:
				del self._callback
				self._callback = NULL
				self._py_callback = None
				raise

	def __dealloc__(self):
		if self._callback != NULL:
			del self._callback
			self._callback = NULL

		if self._parallelCallback != NULL:
			del self._parallelCallback
			self._parallelCallback = NULL

	def setStoreMatches(self, bool_t storeMatches):
		"""
		setStoreMatches(storeMatches)

		Sets whether matches are stored. Pass False to only count them; getMatches() then raises,
		while numberOfMatches() and hasMatch() keep working.

		Parameters
		----------
		storeMatches : bool
			Whether to store matches for getMatches().
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setStoreMatches(storeMatches)

	def getMatches(self):
		"""
		getMatches()

		Returns all matches. Raises a RuntimeError if the matches were not stored, because a
		callback was set or setStoreMatches(False) was called.

		Returns
		-------
		list(list(int))
			The matches, each indexed by pattern node.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).getMatches()

	def numberOfMatches(self):
		"""
		numberOfMatches()

		Returns the number of matches found, whether or not they were stored.

		Returns
		-------
		int
			The number of matches.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).numberOfMatches()

	def hasMatch(self):
		"""
		hasMatch()

		Returns whether at least one match was found.

		Returns
		-------
		bool
			True if at least one match was found.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).hasMatch()


cdef extern from "<networkit/isomorphism/VF2.hpp>" namespace "NetworKit":
	cdef cppclass _VF2 "NetworKit::VF2" (_SubgraphIsomorphism):
		_VF2(const _Graph &pattern, const _Graph &target, _Semantics semantics, count maxMatches) except +

cdef class VF2(SubgraphIsomorphism):
	"""
	VF2(pattern, target, semantics=networkit.isomorphism.Semantics.INDUCED,
		maxMatches=0)

	Finds every occurrence of a pattern graph inside a target graph using the VF2 algorithm.

	Parameters
	----------
	pattern : networkit.Graph
		The graph to look for. Must not contain self-loops.
	target : networkit.Graph
		The graph to look in. Must agree with pattern on directedness.
	semantics : networkit.isomorphism.Semantics, optional
		Whether matches must be induced. Default: networkit.isomorphism.Semantics.INDUCED
	maxMatches : int, optional
		Stop after this many matches; 0 means no limit. Default: 0
	"""

	cdef Graph _pattern
	cdef Graph _target

	def __cinit__(self, Graph pattern, Graph target, semantics=_Semantics.INDUCED, maxMatches=0):
		self._pattern = pattern
		self._target = target
		self._this = new _VF2(pattern._this, target._this, semantics, maxMatches)

cdef extern from "<networkit/isomorphism/RI.hpp>" namespace "NetworKit":
	cdef enum class _RIVariant "NetworKit::RI::Variant":
		RI,
		RI_DS

	cdef cppclass _RI "NetworKit::RI" (_SubgraphIsomorphism):
		_RI(const _Graph &pattern, const _Graph &target, _RIVariant variant, _Semantics semantics, count maxMatches) except +

class Variant(object):
	"""
	Variant of RI and ParallelRI.

	- RI: Plain RI.
	- RI_DS: RI-DS-SI-FC, which computes candidate domains with forward checking before the search.
	"""
	RI = _RIVariant.RI
	RI_DS = _RIVariant.RI_DS

cdef class RI(SubgraphIsomorphism):
	"""
	RI(pattern, target, variant=networkit.isomorphism.Variant.RI,
	   semantics=networkit.isomorphism.Semantics.INDUCED, maxMatches=0)

	Finds every occurrence of a pattern graph inside a target graph using the RI algorithm.

	Bonnici, V., Giugno, R., Pulvirenti, A., Shasha, D., & Ferro, A. (2013).
	A subgraph isomorphism algorithm and its application to biochemical data.
	BMC Bioinformatics, 14(Suppl 7), S13.

	Parameters
	----------
	pattern : networkit.Graph
		The graph to look for. Must not contain self-loops.
	target : networkit.Graph
		The graph to look in. Must agree with pattern on directedness.
	variant : networkit.isomorphism.Variant, optional
		Plain RI or RI-DS. Default: networkit.isomorphism.Variant.RI
	semantics : networkit.isomorphism.Semantics, optional
		Whether matches must be induced. Default: networkit.isomorphism.Semantics.INDUCED
	maxMatches : int, optional
		Stop after this many matches; 0 means no limit. Default: 0
	"""

	cdef Graph _pattern
	cdef Graph _target

	def __cinit__(self, Graph pattern, Graph target, variant=_RIVariant.RI, semantics=_Semantics.INDUCED, maxMatches=0):
		self._pattern = pattern
		self._target = target
		self._this = new _RI(pattern._this, target._this, variant, semantics, maxMatches)

cdef extern from "<networkit/isomorphism/ParallelRI.hpp>":
	cdef cppclass _ParallelRI "NetworKit::ParallelRI" (_SubgraphIsomorphism):
		_ParallelRI(const _Graph &pattern, const _Graph &target, _RIVariant variant, _Semantics semantics, count maxMatches) except +

cdef class ParallelRI(SubgraphIsomorphism):
	"""
	ParallelRI(pattern, target, variant=networkit.isomorphism.Variant.RI,
			   semantics=networkit.isomorphism.Semantics.INDUCED, maxMatches=0)

	Parallel version of RI. Finds the same matches as RI, but in no particular order. The number
	of workers is the global thread count, see networkit.setNumberOfThreads().

	Kimmig, R., Meyerhenke, H., & Strash, D. (2017).
	Shared Memory Parallel Subgraph Enumeration.
	IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW).

	Parameters
	----------
	pattern : networkit.Graph
		The graph to look for. Must not contain self-loops.
	target : networkit.Graph
		The graph to look in. Must agree with pattern on directedness.
	variant : networkit.isomorphism.Variant, optional
		Plain RI or RI-DS. Default: networkit.isomorphism.Variant.RI
	semantics : networkit.isomorphism.Semantics, optional
		Whether matches must be induced. Default: networkit.isomorphism.Semantics.INDUCED
	maxMatches : int, optional
		Stop after this many matches; 0 means no limit. Default: 0
	"""

	cdef Graph _pattern
	cdef Graph _target

	def __cinit__(self, Graph pattern, Graph target, variant=_RIVariant.RI, semantics=_Semantics.INDUCED, maxMatches=0):
		self._pattern = pattern
		self._target = target
		self._this = new _ParallelRI(
			pattern._this, target._this, variant, semantics, maxMatches
		)
