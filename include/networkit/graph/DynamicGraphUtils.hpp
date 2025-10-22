#ifndef NETWORKIT_GRAPH_DYNAMIC_GRAPH_UTILS_HPP_
#define NETWORKIT_GRAPH_DYNAMIC_GRAPH_UTILS_HPP_

#include <algorithm>
#include <tuple>

#include <networkit/Globals.hpp>

namespace NetworKit {

struct Unsafe {};
static constexpr Unsafe unsafe{};

template <class NodeType>
struct Edge {
    NodeType u, v;

    Edge() : u(none), v(none) {}

    Edge(NodeType _u, NodeType _v, bool sorted = false) {
        if (sorted) {
            std::tie(u, v) = std::minmax(_u, _v);
        } else {
            u = _u;
            v = _v;
        }
    }

private:
    static constexpr NodeType none = std::numeric_limits<NodeType>::max();
};

/**
 * A weighted edge used for the graph constructor with
 * initializer list syntax.
 */
template <class NodeType, class EdgeWeightType>
struct WeightedEdge : Edge<NodeType> {
    EdgeWeightType weight;

    // Needed by cython
    WeightedEdge() : Edge<NodeType>(), weight(std::numeric_limits<EdgeWeightType>::max()) {}

    WeightedEdge(NodeType u, NodeType v, EdgeWeightType w) : Edge<NodeType>(u, v), weight(w) {}
};

template <class NodeType>
inline bool operator==(const Edge<NodeType> &e1, const Edge<NodeType> &e2) {
    return e1.u == e2.u && e1.v == e2.v;
}

template <class NodeType, class EdgeWeightType>
struct WeightedEdgeWithId : WeightedEdge<NodeType, EdgeWeightType> {
    edgeid eid;

    WeightedEdgeWithId(NodeType u, NodeType v, EdgeWeightType w, edgeid eid)
        : WeightedEdge<NodeType, EdgeWeightType>(u, v, w), eid(eid) {}
};

template <class NodeType, class EdgeWeightType>
inline bool operator<(const WeightedEdge<NodeType, EdgeWeightType> &e1,
                      const WeightedEdge<NodeType, EdgeWeightType> &e2) {
    return e1.weight < e2.weight;
}

// We expose to Cython the default 64-bit Edge and WeightedEdge.
using _CythonEdge = Edge<node>;
using _CythonWeightedEdge = WeightedEdge<node, edgeweight>;
} // namespace NetworKit

namespace std {
template <class NodeType>
struct hash<NetworKit::Edge<NodeType>> {
    size_t operator()(const NetworKit::Edge<NodeType> &edge) const {
        return hash_node(edge.u) ^ hash_node(edge.v);
    }
    hash<NodeType> hash_node;
};
} // namespace std

#endif // NETWORKIT_GRAPH_DYNAMIC_GRAPH_UTILS_HPP_
