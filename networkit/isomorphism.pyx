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
	INDUCED = _Semantics.INDUCED
	MONOMORPHISM = _Semantics.MONOMORPHISM

cdef class SubgraphIsomorphism(Algorithm):
	"""Abstract base class for subgraph isomorphism algorithms."""

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
		Return the number of workers used by the algorithm.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).numberOfWorkers()

	def setNodeLabels(self, vector[index] patternNodeLabels, vector[index] targetNodeLabels):
		"""
		setNodeLabels(patternNodeLabels, targetNodeLabels)

		Set labels for pattern and target nodes. Labels are indexed by node id.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setNodeLabels(patternNodeLabels, targetNodeLabels)

	def setEdgeLabels(self, vector[index] patternEdgeLabels,
					  vector[index] targetEdgeLabels):
		"""
		setEdgeLabels(patternEdgeLabels, targetEdgeLabels)

		Set labels for pattern and target edges. Labels are indexed by edge id.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setEdgeLabels(patternEdgeLabels, targetEdgeLabels)

	def setCallback(self, object callback, bool_t parallel=False):
		"""
		setCallback(callback, parallel=False)

		Set a Python callback that is called for every match.

		With ``parallel=False``:
			callback(match)

		With ``parallel=True``:
			callback(workerId, match)

		If the callback raises, the search stops and ``run()`` raises a
		``RuntimeError`` that carries the original message.
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

		Choose whether found matches are kept for ``getMatches()``.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).setStoreMatches(storeMatches)

	def getMatches(self):
		"""
		getMatches()

		Return all stored matches. Each match is indexed by pattern node.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).getMatches()

	def numberOfMatches(self):
		"""
		Return the number of matches found.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return (<_SubgraphIsomorphism*>(self._this)).numberOfMatches()

	def hasMatch(self):
		"""
		Return whether at least one match was found.
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

	Finds every occurrence of a pattern graph inside a target graph using VF2.
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
	RI = _RIVariant.RI
	RI_DS = _RIVariant.RI_DS

cdef class RI(SubgraphIsomorphism):
	"""
	RI(pattern, target, variant=networkit.isomorphism.Variant.RI,
	   semantics=networkit.isomorphism.Semantics.INDUCED, maxMatches=0)

	Finds every occurrence of a pattern graph inside a target graph using RI.
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

	Parallel version of RI.
	"""

	cdef Graph _pattern
	cdef Graph _target

	def __cinit__(self, Graph pattern, Graph target, variant=_RIVariant.RI, semantics=_Semantics.INDUCED, maxMatches=0):
		self._pattern = pattern
		self._target = target
		self._this = new _ParallelRI(
			pattern._this, target._this, variant, semantics, maxMatches
		)