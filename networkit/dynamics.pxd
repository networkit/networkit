from libcpp cimport bool as bool_t
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair

from .structures cimport node, index, count, edgeweight

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
		bool_t operator==(_GraphEvent)
		bool_t operator!=(_GraphEvent)
		bool_t operator<(_GraphEvent, _GraphEvent)
		bool_t operator>(_GraphEvent, _GraphEvent)
		bool_t operator<=(_GraphEvent)
		bool_t operator>=(_GraphEvent)

cdef class GraphEvent:
	cdef _GraphEvent _this
