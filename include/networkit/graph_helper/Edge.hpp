#ifndef NETWORKIT_GRAPH_HELPER_EDGE
#define NETWORKIT_GRAPH_HELPER_EDGE

#include <algorithm>
#include <vector>
#include <networkit/Globals.hpp>

namespace NetworKit {

struct Edge {
    node u, v;

    Edge() : u(none), v(none) {}

    Edge(node _u, node _v, bool sorted = false) {
        u = sorted ? std::min(_u, _v) : _u;
        v = sorted ? std::max(_u, _v) : _v;
    }
};

/**
 * A weighted edge used for the graph constructor with
 * initializer list syntax.
 */
struct WeightedEdge : Edge {
    edgeweight weight;

    // Needed by cython
    WeightedEdge() : Edge(), weight(std::numeric_limits<edgeweight>::max()) {}

    WeightedEdge(node u, node v, edgeweight w) : Edge(u, v), weight(w) {}
};

struct WeightedEdgeWithId : WeightedEdge {
    edgeid eid;

    WeightedEdgeWithId(node u, node v, edgeweight w, edgeid eid)
        : WeightedEdge(u, v, w), eid(eid) {}
};

inline bool operator==(const Edge &e1, const Edge &e2) {
    return e1.u == e2.u && e1.v == e2.v;
}

inline bool operator<(const WeightedEdge &e1, const WeightedEdge &e2) {
    return e1.weight < e2.weight;
}

struct Unsafe {};
static constexpr Unsafe unsafe{};
} // namespace NetworKit

namespace std {
template <>
struct hash<NetworKit::Edge> {
    size_t operator()(const NetworKit::Edge &e) const { return hash_node(e.u) ^ hash_node(e.v); }

    hash<NetworKit::node> hash_node;
};
} // namespace std

#endif