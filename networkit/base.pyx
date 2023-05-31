# distutils: language=c++


cdef class _CythonParentClass:
	""" Abstract base class for Cython 
	The purpose of this class is to provide a single combined base cdef class for all classes that we use in this project.
	Cython cannot handle multiple inheritance without this single common base class.
	"""
	def __init__(self, *args, **namedargs):
		if type(self) == _CythonParentClass:
			raise RuntimeError("Error, you may not use _CythonParentClass directly, use a sub-class instead")

	def __cinit__(self, *args, **namedargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

cdef class Algorithm(_CythonParentClass):
	""" Abstract base class for algorithms """
	def __init__(self, *args, **namedargs):
		if type(self) == Algorithm:
			raise RuntimeError("Error, you may not use Algorithm directly, use a sub-class instead")

	def __cinit__(self, *args, **namedargs):
		self._this = NULL

	def __dealloc__(self):
		if self._this != NULL:
			del self._this
		self._this = NULL

	def run(self):
		"""
		run()

		Executes the algorithm.

		Returns
		-------
		networkit.base.Algorithm
			self
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		with nogil:
			self._this.run()
		return self

	def hasFinished(self):
		"""
		hasFinished()

		States whether an algorithm has already run.

		Returns
		-------
		bool
			True if Algorithm has finished.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.hasFinished()
