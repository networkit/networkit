#ifndef NETWORKIT_GRAPH_EDGE_UTILS_HPP_
#define NETWORKIT_GRAPH_EDGE_UTILS_HPP_

#include <algorithm>
#include <tuple>

#include <networkit/Globals.hpp>

namespace NetworKit {

struct Unsafe {};
static constexpr Unsafe unsafe{};

template <class NodeT>
struct EdgeT {
    NodeT u, v;

    EdgeT() : u(nullNodeId), v(nullNodeId) {}

    EdgeT(NodeT _u, NodeT _v, bool sorted = false) {
        if (sorted) {
            std::tie(u, v) = std::minmax(_u, _v);
        } else {
            u = _u;
            v = _v;
        }
    }

private:
    static constexpr NodeT nullNodeId = NullNodeId<NodeT>;
};

/**
 * A weighted edge used for the graph constructor with
 * initializer list syntax.
 */
template <class NodeT, class EdgeWeightT>
struct WeightedEdgeT : EdgeT<NodeT> {
    EdgeWeightT weight;

    // Needed by cython
    WeightedEdgeT() : EdgeT<NodeT>(), weight(std::numeric_limits<EdgeWeightT>::max()) {}

    WeightedEdgeT(NodeT u, NodeT v, EdgeWeightT w) : EdgeT<NodeT>(u, v), weight(w) {}
};

template <class NodeT>
inline bool operator==(const EdgeT<NodeT> &e1, const EdgeT<NodeT> &e2) {
    return e1.u == e2.u && e1.v == e2.v;
}

template <class NodeT, class EdgeWeightT>
struct WeightedEdgeWithId : WeightedEdgeT<NodeT, EdgeWeightT> {
    edgeid eid;

    WeightedEdgeWithId(NodeT u, NodeT v, EdgeWeightT w, edgeid eid)
        : WeightedEdgeT<NodeT, EdgeWeightT>(u, v, w), eid(eid) {}
};

template <class NodeT, class EdgeWeightT>
inline bool operator<(const WeightedEdgeT<NodeT, EdgeWeightT> &e1,
                      const WeightedEdgeT<NodeT, EdgeWeightT> &e2) {
    return e1.weight < e2.weight;
}

using Edge = EdgeT<node>;
using WeightedEdge = WeightedEdgeT<node, edgeweight>;

// We expose to Cython the default 64-bit Edge and WeightedEdge.
using _CythonEdge = Edge;
using _CythonWeightedEdge = WeightedEdge;
} // namespace NetworKit

namespace std {
template <class NodeT>
struct hash<NetworKit::EdgeT<NodeT>> {
    size_t operator()(const NetworKit::EdgeT<NodeT> &edge) const {
        return hash_node(edge.u) ^ hash_node(edge.v);
    }
    hash<NodeT> hash_node;
};
} // namespace std

#endif // NETWORKIT_GRAPH_EDGE_UTILS_HPP_
