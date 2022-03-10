# distutils: language=c++

cdef class Algorithm:
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
