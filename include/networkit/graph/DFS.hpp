
#ifndef NETWORKIT_GRAPH_DFS_HPP_
#define NETWORKIT_GRAPH_DFS_HPP_

#include <cstddef>
#include <stack>
#include <stdexcept>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/TraversalImpl.hpp>

namespace NetworKit {

namespace Traversal {

/**
 * Iterate over nodes in depth-first search order starting from the given source node.
 *
 * @param G The input graph.
 * @param source The source node.
 * @param handle Takes a node as input parameter.
 * @throws std::runtime_error If the source node is not in the graph.
 */
template <typename L, class GraphT>
void DFSfrom(const GraphT &G, typename GraphT::NodeT source, L handle) {
    using NodeT = typename GraphT::NodeT;

    if (!Impl::isValidNode(G, source))
        throw std::runtime_error("Error: source node not in the graph.");
    std::vector<Visit> marked(static_cast<std::size_t>(G.upperNodeIdBound()), Visit::Unvisited);
    std::stack<NodeT> s;
    s.push(source); // enqueue root
    marked[Impl::nodeIndex(source)] = Visit::Visited;
    do {
        const auto u = s.top();
        s.pop();
        // apply function
        handle(u);
        G.forNeighborsOf(u, [&](NodeT v) {
            const std::size_t vIndex = Impl::nodeIndex(v);
            if (marked[vIndex] == Visit::Unvisited) {
                s.push(v);
                marked[vIndex] = Visit::Visited;
            }
        });
    } while (!s.empty());
}

/**
 * Iterate over edges in depth-first search order starting from the given source node.
 *
 * @param G The input graph.
 * @param source The source node.
 * @param handle Takes the source node, target node, weight and ID of each discovery edge.
 * @throws std::runtime_error If the source node is not in the graph.
 */
template <typename L, class GraphT>
void DFSEdgesFrom(const GraphT &G, typename GraphT::NodeT source, L handle) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    if (!Impl::isValidNode(G, source))
        throw std::runtime_error("Error: source node not in the graph.");
    std::vector<Visit> marked(static_cast<std::size_t>(G.upperNodeIdBound()), Visit::Unvisited);
    std::stack<NodeT> s;
    s.push(source); // enqueue root
    marked[Impl::nodeIndex(source)] = Visit::Visited;
    do {
        const auto u = s.top();
        s.pop();
        // apply function
        G.forNeighborsOf(u, [&](NodeT, NodeT v, EdgeWeightT w, edgeid eid) {
            const std::size_t vIndex = Impl::nodeIndex(v);
            if (marked[vIndex] == Visit::Unvisited) {
                handle(u, v, w, eid);
                s.push(v);
                marked[vIndex] = Visit::Visited;
            }
        });
    } while (!s.empty());
}

} // namespace Traversal

} // namespace NetworKit
#endif // NETWORKIT_GRAPH_DFS_HPP_
