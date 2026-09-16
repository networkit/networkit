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

		Set labels for pattern and target nodes. Labels are indexed by node id.

		Parameters
		----------
		patternNodeLabels : list(int)
			List with pattern node labels.
		targetNodeLabels : list(int)
			List with target node labels.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setNodeLabels(patternNodeLabels, targetNodeLabels)

	def setEdgeLabels(self, vector[index] patternEdgeLabels,
					  vector[index] targetEdgeLabels):
		"""
		setEdgeLabels(patternEdgeLabels, targetEdgeLabels)

		Set labels for pattern and target edges. Labels are indexed by edge id.

		Parameters
		----------
		patternEdgeLabels : list(int)
			List with pattern edge labels.
		targetEdgeLabels : list(int)
			List with target edge labels.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setEdgeLabels(patternEdgeLabels, targetEdgeLabels)

	def setSequentialCallback(self, object callback):
		"""
		setSequentialCallback(callback)

		Set a Python callback. Every match will be handed to this callback as it is found, rather than collecting them. If the callback raises, the search stops and ``run()`` raises a ``RuntimeError`` that carries the original message.

		Parameters
		----------
		callback : callable
			Called once per match. Must accept (match).
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")

		if not callable(callback):
			raise TypeError("Callback must be callable")

		if isinstance(self, ParallelRI):
			raise RuntimeError("Error, must use setParallelCallback(callback) for parallel algorithms")

		try:
			signature = inspect.signature(callback)
			signature.bind(None)
		except (TypeError, ValueError):
				raise TypeError("callback must accept (match)") from None

		# Remove a previously registered wrapper.
		if self._callback != NULL:
			del self._callback
			self._callback = NULL

		# Keep the Python callback alive.
		self._py_callback = callback

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

	def setParallelCallback(self, object callback):
		"""
		setParallelCallback(callback)

		Like setSequentialCallback(callback), but this one sets a callback that also receives the worker id.

		Parameters
		----------
		callback : callable
			Called once per match. Must accept (workerId, match).
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")

		if not callable(callback):
			raise TypeError("Callback must be callable")

		if (isinstance(self, ParallelRI) == False):
			raise RuntimeError ("Error, must use setSequentialCallback(callback) for sequential algorithms")

		try:
			signature = inspect.signature(callback)
			signature.bind(None, None)
		except (TypeError, ValueError):
			raise TypeError("callback must accept (workerId, match)") from None

		# Remove a previously registered wrapper.
		if self._parallelCallback != NULL:
			del self._parallelCallback
			self._parallelCallback = NULL

		# Keep the Python callback alive.
		self._py_callback = callback

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

		Choose whether found matches are kept for getMatches().

		Parameters
		----------
		storeMatches : bool
			Whether to keep matches for getMatches(). Default: True
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setStoreMatches(storeMatches)

	def getMatches(self):
		"""
		getMatches()

		Returns all matches found, each one a vector indexed by pattern node. Throws an error if the matches were never stored, which happens when a callback was set or @ref setStoreMatches(false) was called, and if @ref run() has not been called yet.

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
		Returns how many matches were found. Works regardless of whether they were stored. If a match limit was specified, the returned value is capped at that limit.

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

		Return whether at least one match was found.
	
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

	Finds every occurrence of a pattern graph inside a target graph, using the VF2 algorithm.

	Parameters
	----------
	pattern : networkit.Graph
		The pattern graph. Must not contain self-loops.
	target : networkit.Graph
		The target graph. Must agree with pattern on directedness.
	semantics : 
		Whether matches must be induced. Default: networkit.isomorphism.Semantics.INDUCED
	maxMatches : int
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

	Finds every occurrence of a pattern graph inside a target graph using RI.

	Parameters
	----------
	pattern : networkit.Graph
		The pattern graph. Must not contain self-loops.
	target : networkit.Graph
		The target graph. Must agree with pattern on directedness.
	variant : networkit.isomorphism.Variant.RI
		Plain RI or RI-DS.
	semantics : 
		Whether matches must be induced. Default: networkit.isomorphism.Semantics.INDUCED
	maxMatches : int
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

	Parallel version of RI.
	
	Parameters
	----------
	pattern : networkit.Graph
		The pattern graph. Must not contain self-loops.
	target : networkit.Graph
		The target graph. Must agree with pattern on directedness.
	variant : networkit.isomorphism.Variant.RI
		Plain RI or RI-DS.
	semantics : 
		Whether matches must be induced. Default: networkit.isomorphism.Semantics.INDUCED
	maxMatches : int
		Stop after this many matches; 0 means no limit. Default: 0
	"""

	def __cinit__(self, Graph pattern, Graph target, variant=_RIVariant.RI, semantics=_Semantics.INDUCED, maxMatches=0):
		self._pattern = pattern
		self._target = target
		self._this = new _ParallelRI(
			pattern._this, target._this, variant, semantics, maxMatches
		)
