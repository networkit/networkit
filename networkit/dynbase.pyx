# distutils: language=c++

from libcpp.cast cimport dynamic_cast

from .dynamics cimport _GraphEvent, GraphEvent
from .base cimport _CythonParentClass
from typing import List

ctypedef _DynAlgorithm* _DynAlgorithmPtr

cdef class _CythonSubclassDynAlgorithm(_CythonParentClass):
	""" Abstract base class for dynamic algorithms
	Cython does not allow direct inheritance from two extension types (cdef classes).
	Instead, inherit from the DynAlgorithm class! 
	"""
	def __init__(self, *args, **namedargs):
		if type(self) == _CythonSubclassDynAlgorithm:
			raise RuntimeError("Error, you may not use _CythonSubclassDynAlgorithm directly, use a sub-class instead")

	def __cinit__(self, *args, **namedargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

	def update(self, GraphEvent ev):
		"""
		update(ev)

		The generic update method for updating data structure after an update.

		Parameters
		----------
		ev : networkit.dynamics.GraphEvent
			A graph event.

		Returns
		-------
		networkit.dynbase.DynAlgorithm
			self
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		cdef _GraphEvent _ev = _GraphEvent(ev.type, ev.u, ev.v, ev.w)
		with nogil:
			dynamic_cast[_DynAlgorithmPtr](self._this).update(_ev)
		return self
	
	def updateBatch(self, batch: List[GraphEvent]):
		""" 
		updateBatch(batch)

		The generic update method for updating data structure after a batch of updates.

		Parameters
		----------
		batch : list(networkit.dynamics.GraphEvent)
			List of graph events.

		Returns
		-------
		networkit.dynbase.DynAlgorithm
			self
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		cdef vector[_GraphEvent] _batch
		for ev in batch:
			_batch.push_back(_GraphEvent(ev.type, ev.u, ev.v, ev.w))
		with nogil:
			dynamic_cast[_DynAlgorithmPtr](self._this).updateBatch(_batch)
		return self

class DynAlgorithm(_CythonSubclassDynAlgorithm):
	__slots__ = ()
	""" Abstract base class for dynamic algorithms """
	def __init__(self, *args, **namedargs):
		if type(self) == DynAlgorithm:
			raise RuntimeError("Error, you may not use DynAlgorithm directly, use a sub-class instead")
