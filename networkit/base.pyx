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
		Executes the algorithm.

		Returns:
		--------
		Algorithm
			self
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		with nogil:
			self._this.run()
		return self

	def hasFinished(self):
		"""
		States whether an algorithm has already run.

		Returns:
		--------
		bool
			True if Algorithm has finished.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.hasFinished()

	def toString(self):
		""" Get string representation.

		Returns:
		--------
		string
			String representation of algorithm and parameters.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.toString().decode("utf-8")


	def isParallel(self):
		"""
		Returns:
		--------
		bool
			True if algorithm can run multi-threaded.
		"""
		if self._this == NULL:
			raise RuntimeError("Error, object not properly initialized")
		return self._this.isParallel()
