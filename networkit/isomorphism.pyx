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

import inspect

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
			except BaseException as e:
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
			except BaseException as e:
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
	non-edges. The search finds an occurrence once per automorphism of the pattern, so a triangle
	occurs six times.

	VF2 is the reference implementation, and RI is usually faster on sparse targets. ParallelRI runs
	RI on several threads and reports the matches in no fixed order.

	Pattern and target must both be directed or both be undirected, and the pattern must not contain
	self-loops. The search ignores edge weights and target self-loops, and it collapses parallel
	edges.
	"""

	cdef Graph _pattern
	cdef Graph _target
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
		networkit.none matches any label, and two empty lists remove the labels. Call this before
		run().

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
		networkit.none matches any label, and two empty lists remove the labels. Both graphs need
		edge ids, see networkit.Graph.indexEdges(). run() raises an error for parallel edges with
		different labels. Call this before run().

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

	def setCallback(self, object callback):
		"""
		setCallback(callback)

		Passes every match to the callback instead of storing it, so getMatches() raises an error.
		The callback is never called concurrently. ParallelRI makes its
		workers take turns at the callback, so the callback may change shared state without a lock.
		If the callback raises, the search stops and ``run()`` raises a ``RuntimeError`` that
		carries the original message. A later call of setCallback() or setParallelCallback()
		replaces the callback.

		Parameters
		----------
		callback : callable
			Called once per match. Must accept (match).
		"""
		cdef MatchCallbackWrapper* wrapper

		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")

		if not callable(callback):
			raise TypeError("Callback must be callable")

		try:
			signature = inspect.signature(callback)
			signature.bind(None)
		except (TypeError, ValueError):
			raise TypeError("callback must accept (match)") from None

		# setCallback() stores its own copy of the wrapper, so the wrapper is freed right away.
		wrapper = new MatchCallbackWrapper(callback)
		try:
			(<_SubgraphIsomorphism*>(self._this)).setCallback(dereference(wrapper))
		finally:
			del wrapper

		# The copy holds no reference to the Python callback, so keep the callback alive here.
		self._py_callback = callback

	def setParallelCallback(self, object callback):
		"""
		setParallelCallback(callback)

		Like setCallback(callback), but the callback also receives the id of the worker that found
		the match. ParallelRI may call the callback from several workers, and these calls can
		overlap, so the callback must be thread-safe. The worker id is smaller than
		numberOfWorkers(), so results that are kept per worker id need no lock. ParallelRI uses the
		global thread count, so call numberOfWorkers() after networkit.setNumberOfThreads(). The
		sequential algorithms always pass worker id 0.

		Parameters
		----------
		callback : callable
			Called once per match. Must accept (workerId, match).
		"""
		cdef ParallelMatchCallbackWrapper* wrapper

		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")

		if not callable(callback):
			raise TypeError("Callback must be callable")

		try:
			signature = inspect.signature(callback)
			signature.bind(None, None)
		except (TypeError, ValueError):
			raise TypeError("callback must accept (workerId, match)") from None

		# setCallback() stores its own copy of the wrapper, so the wrapper is freed right away.
		wrapper = new ParallelMatchCallbackWrapper(callback)
		try:
			(<_SubgraphIsomorphism*>(self._this)).setCallback(dereference(wrapper))
		finally:
			del wrapper

		# The copy holds no reference to the Python callback, so keep the callback alive here.
		self._py_callback = callback

	def setStoreMatches(self, bool_t storeMatches):
		"""
		setStoreMatches(storeMatches)

		Sets whether matches are stored. Pass False to only count them, so getMatches() raises an
		error. Matches are stored by default. Call this before run().

		Parameters
		----------
		storeMatches : bool
			Whether to keep matches for getMatches().
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setStoreMatches(storeMatches)

	def getMatches(self):
		"""
		getMatches()

		Returns all matches found, each one a vector indexed by pattern node. Raises an error if
		run() has not been called yet or if the matches were not stored. Matches are not stored when
		a callback was set or setStoreMatches(False) was called.

		Returns
		-------
		list(list(int))
			A list of matches, each being represented as a list of target node ids.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).getMatches()

	def numberOfMatches(self):
		"""
		numberOfMatches()

		Returns how many matches were found, stored or not. The returned value is at most
		maxMatches. ParallelRI with a callback is the exception, since its workers may deliver a few
		matches beyond maxMatches before they all stop. These matches are counted too.

		Returns
		-------
		int
			The number of matches reported.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).numberOfMatches()

	def hasMatch(self):
		"""
		hasMatch()

		Returns whether a match was found. Pass maxMatches=1 to stop at the first one.

		Returns
		-------
		bool
			True if there was at least one match.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).hasMatch()

cdef extern from "<networkit/isomorphism/VF2.hpp>" namespace "NetworKit":
	cdef cppclass _VF2 "NetworKit::VF2" (_SubgraphIsomorphism):
		_VF2(const _Graph &pattern, const _Graph &target, _Semantics semantics, count maxMatches) except +

cdef class VF2(SubgraphIsomorphism):
	"""
	VF2(pattern, target, semantics=networkit.isomorphism.Semantics.INDUCED, maxMatches=0)

	Finds every occurrence of a pattern graph inside a target graph using the VF2 algorithm.

	VF2 extends a partial mapping one node pair at a time, depth first, and draws the candidate
	pairs from the terminal sets of unmapped nodes adjacent to mapped ones. VF2 serves as the
	reference implementation, but RI is usually faster on sparse targets.

	Parameters
	----------
	pattern : networkit.Graph
		The pattern graph. Must not contain self-loops.
	target : networkit.Graph
		The target graph. Must agree with pattern on directedness.
	semantics : networkit.isomorphism.Semantics, optional
		Whether matches must be induced. Default: networkit.isomorphism.Semantics.INDUCED
	maxMatches : int, optional
		Stop after this many matches; 0 means no limit. Default: 0
	"""

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

	RI fixes the order of the pattern nodes once and then backtracks along this order. The variant
	networkit.isomorphism.Variant.RI_DS first computes a domain of candidate target nodes for every
	pattern node. It pays off for disconnected patterns and selective node labels, but otherwise
	usually costs more than it saves.

	Parameters
	----------
	pattern : networkit.Graph
		The pattern graph. Must not contain self-loops.
	target : networkit.Graph
		The target graph. Must agree with pattern on directedness.
	variant : networkit.isomorphism.Variant, optional
		Plain RI or RI-DS. Default: networkit.isomorphism.Variant.RI
	semantics : networkit.isomorphism.Semantics, optional
		Whether matches must be induced. Default: networkit.isomorphism.Semantics.INDUCED
	maxMatches : int, optional
		Stop after this many matches; 0 means no limit. Default: 0
	"""

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

	Parallel version of RI. It finds the same matches as RI, but their order may differ from run to
	run. The number of workers is the global thread count, see networkit.setNumberOfThreads(). Idle
	workers steal batches of partial mappings from busy ones. A callback set with setCallback()
	makes the workers take turns, and a callback set with setParallelCallback() avoids this.

	Parameters
	----------
	pattern : networkit.Graph
		The pattern graph. Must not contain self-loops.
	target : networkit.Graph
		The target graph. Must agree with pattern on directedness.
	variant : networkit.isomorphism.Variant, optional
		Plain RI or RI-DS. Default: networkit.isomorphism.Variant.RI
	semantics : networkit.isomorphism.Semantics, optional
		Whether matches must be induced. Default: networkit.isomorphism.Semantics.INDUCED
	maxMatches : int, optional
		Stop after this many matches; 0 means no limit. Default: 0
	"""

	def __cinit__(self, Graph pattern, Graph target, variant=_RIVariant.RI, semantics=_Semantics.INDUCED, maxMatches=0):
		self._pattern = pattern
		self._target = target
		self._this = new _ParallelRI(
			pattern._this, target._this, variant, semantics, maxMatches
		)
