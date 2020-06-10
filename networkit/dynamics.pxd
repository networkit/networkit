from libc.stdint cimport uint64_t
from libcpp cimport bool as bool_t
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair

ctypedef uint64_t count
ctypedef uint64_t index
ctypedef index node
ctypedef double edgeweight

cdef extern from "<networkit/dynamics/GraphEvent.hpp>" namespace "NetworKit::GraphEvent::Type":

	cdef enum _GraphEventType "NetworKit::GraphEvent::Type":
		NODE_ADDITION,
		NODE_REMOVAL,
		NODE_RESTORATION,
		EDGE_ADDITION,
		EDGE_REMOVAL,
		EDGE_WEIGHT_UPDATE,
		EDGE_WEIGHT_INCREMENT,
		TIME_STEP

cdef extern from "<networkit/dynamics/GraphEvent.hpp>":

	cdef cppclass _GraphEvent "NetworKit::GraphEvent":
		node u, v
		edgeweight w
		_GraphEventType type
		_GraphEvent() except +
		_GraphEvent(_GraphEventType type, node u, node v, edgeweight w) except +
		string toString() except +

cdef class GraphEvent:
	cdef _GraphEvent _this

cdef extern from "<networkit/dynamics/GraphEvent.hpp>" namespace "NetworKit::GraphEvent":

	bool_t _GraphEvent_equal "NetworKit::GraphEvent::equal"(_GraphEvent a, _GraphEvent b) except +
	bool_t _GraphEvent_compare "NetworKit::GraphEvent::compare"(_GraphEvent a, _GraphEvent b) except +

