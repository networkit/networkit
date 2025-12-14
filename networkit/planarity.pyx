# distutils: language = c++

from libcpp cimport bool as cbool

from .graph cimport _Graph, Graph
from .base cimport _Algorithm, Algorithm


cdef extern from "<networkit/planarity/LeftRightPlanarityCheck.hpp>":
    cdef cppclass _LeftRightPlanarityCheck "NetworKit::LeftRightPlanarityCheck"(_Algorithm):
        _LeftRightPlanarityCheck(const _Graph& G) except +
        void run() except +
        cbool isPlanar() except +

cdef class LeftRightPlanarityCheck(Algorithm):
    """
    Left-Right Planarity Test.

    Parameters
    ----------
    G : networkit.Graph
        The input graph.
    """
    cdef Graph _graph
    cdef _LeftRightPlanarityCheck* _this

    def __cinit__(self, Graph graph not None):
        self._graph = graph
        self._this = new _LeftRightPlanarityCheck(graph._this)

    def isPlanar(self):
        """
        Returns True iff the graph is planar.

        """
        return bool(self._this.isPlanar())
