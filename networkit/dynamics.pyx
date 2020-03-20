cdef class GraphEvent:
	NODE_ADDITION = 0
	NODE_REMOVAL = 1
	NODE_RESTORATION = 2
	EDGE_ADDITION = 3
	EDGE_REMOVAL = 4
	EDGE_WEIGHT_UPDATE = 5
	EDGE_WEIGHT_INCREMENT = 6
	TIME_STEP = 7

	property type:
		def __get__(self):
			return self._this.type
		def __set__(self, t):
			self._this.type = t

	property u:
		def __get__(self):
			return self._this.u
		def __set__(self, u):
			self._this.u = u

	property v:
		def __get__(self):
			return self._this.v
		def __set__(self, v):
			self._this.v = v

	property w:
		def __get__(self):
			return self._this.w
		def __set__(self, w):
			self._this.w = w

	def __cinit__(self, _GraphEventType type, node u, node v, edgeweight w):
		self._this = _GraphEvent(type, u, v, w)

	def toString(self):
		return self._this.toString().decode("utf-8")

	def __repr__(self):
		return self.toString()

	def __eq__(self, GraphEvent other not None):
		return _GraphEvent_equal(self._this, other._this)
