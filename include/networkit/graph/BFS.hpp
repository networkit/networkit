
#ifndef NETWORKIT_GRAPH_BFS_HPP_
#define NETWORKIT_GRAPH_BFS_HPP_

#include <array>
#include <cassert>
#include <queue>
#include <type_traits>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

namespace Traversal {

namespace Impl {

/**
 * Calls the given BFS handle with distance parameter
 */
template <class F, class NodeT>
auto callBFSHandle(F &f, NodeT u, count dist) -> decltype(f(u, dist)) {
    return f(u, dist);
}

/**
 * Calls the given BFS handle without distance parameter
 */
template <class F, class NodeT>
auto callBFSHandle(F &f, NodeT u, count) -> decltype(f(u)) {
    return f(u);
}

template <typename NodeT>
index nodeIndex(NodeT u) {
    if constexpr (std::is_signed_v<NodeT>) {
        assert(u >= 0);
    }
    return static_cast<index>(u);
}

} // namespace Impl

/**
 * Iterate over nodes in breadth-first search order starting from the unique nodes within the
 * given range. All start nodes are visited at distance 0. Duplicate start nodes are ignored.
 *
 * @param G The input graph.
 * @param first The first element of the start node range.
 * @param last The end of the start node range.
 * @param handle Takes a node, or a node and its distance from the nearest start node.
 */
template <class GraphT, class InputIt, typename L>
void BFSfrom(const GraphT &G, InputIt first, InputIt last, L handle) {
    using NodeT = typename GraphT::NodeT;

    std::vector<bool> marked(static_cast<index>(G.upperNodeIdBound()));
    std::queue<NodeT> q, qNext;
    count dist = 0;
    // enqueue start nodes, do not enqueue duplicates
    for (; first != last; ++first) {
        const index uIndex = Impl::nodeIndex(*first);
        if (!marked[uIndex]) {
            q.push(*first);
            marked[uIndex] = true;
        }
    }
    while (!q.empty()) {
        const auto u = q.front();
        q.pop();
        // apply function
        Impl::callBFSHandle(handle, u, dist);
        G.forNeighborsOf(u, [&](NodeT v) {
            if (!marked[Impl::nodeIndex(v)]) {
                qNext.push(v);
                marked[Impl::nodeIndex(v)] = true;
            }
        });
        if (q.empty() && !qNext.empty()) {
            q.swap(qNext);
            ++dist;
        }
    }
}

/**
 * Iterate over nodes in breadth-first search order starting from the given source node.
 *
 * @param G The input graph.
 * @param source The source node.
 * @param handle Takes a node as input parameter.
 */
template <class GraphT, typename L>
void BFSfrom(const GraphT &G, typename GraphT::NodeT source, L handle) {
    using NodeT = typename GraphT::NodeT;

    std::array<NodeT, 1> startNodes{{source}};
    BFSfrom(G, startNodes.begin(), startNodes.end(), handle);
}

/**
 * Iterate over edges in breadth-first search order starting from the given source node.
 *
 * @param G The input graph.
 * @param source The source node.
 * @param handle Takes a node as input parameter.
 */
template <class GraphT, typename L>
void BFSEdgesFrom(const GraphT &G, typename GraphT::NodeT source, L handle) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    std::vector<bool> marked(static_cast<index>(G.upperNodeIdBound()));
    std::queue<NodeT> q;
    q.push(source); // enqueue root
    marked[Impl::nodeIndex(source)] = true;
    do {
        const auto u = q.front();
        q.pop();
        // apply function
        G.forNeighborsOf(u, [&](NodeT, NodeT v, EdgeWeightT w, edgeid eid) {
            if (!marked[Impl::nodeIndex(v)]) {
                handle(u, v, w, eid);
                q.push(v);
                marked[Impl::nodeIndex(v)] = true;
            }
        });
    } while (!q.empty());
}

} // namespace Traversal

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_BFS_HPP_
