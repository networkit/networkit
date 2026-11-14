
#ifndef NETWORKIT_GRAPH_DIJKSTRA_HPP_
#define NETWORKIT_GRAPH_DIJKSTRA_HPP_

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/TraversalImpl.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {

namespace Traversal {

namespace Impl {
/**
 * Calls the given Dijkstra handle with distance parameter
 */
template <class F, class NodeT, class EdgeWeight>
auto callDijkstraHandle(F &f, NodeT u, EdgeWeight dist) -> decltype(f(u, dist)) {
    return f(u, dist);
}

/**
 * Calls the given Dijkstra handle without distance parameter
 */
template <class F, class NodeT, class EdgeWeight>
auto callDijkstraHandle(F &f, NodeT u, EdgeWeight) -> decltype(f(u)) {
    return f(u);
}
} // namespace Impl

/**
 * Iterate over nodes with Dijkstra starting from the nodes within the given range.
 * Duplicate source nodes are ignored.
 *
 * @param G The input graph.
 * @param first The beginning of the range of source nodes.
 * @param last The end of the range.
 * @param handle Takes a node and optionally its distance from the nodes in the range.
 * @throws std::runtime_error If a source node is not in the graph.
 */
template <class InputIt, typename Handle, class GraphT>
    requires Impl::TraversalHandle<Handle, typename GraphT::NodeT, typename GraphT::EdgeWeightT>
void DijkstraFrom(const GraphT &G, InputIt first, InputIt last, Handle handle) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
    std::vector<EdgeWeightT> distance(G.upperNodeIdBound(),
                                      std::numeric_limits<EdgeWeightT>::max());
    const auto compareDistance = [&distance](std::size_t u, std::size_t v) noexcept -> bool {
        return distance[u] < distance[v];
    };

    tlx::d_ary_addressable_int_heap<std::size_t, 2, decltype(compareDistance)> heap{
        compareDistance};

    std::for_each(first, last, [&](const NodeT u) {
        if (!Impl::isValidNode(G, u))
            throw std::runtime_error("Error: source node not in the graph.");
        // During initialization, zero identifies an already inserted source.
        if (distance[u] == EdgeWeightT{0})
            return;
        distance[u] = EdgeWeightT{0};
        heap.push(Impl::nodeIndex(u));
    });

    while (!heap.empty()) {
        const auto u = static_cast<NodeT>(heap.extract_top());
        Impl::callDijkstraHandle(handle, u, distance[u]);
        G.forNeighborsOf(
            u, [&](NodeT /* source */, const NodeT v, const EdgeWeightT w, edgeid /*id*/) {
                if (distance[v] > distance[u] + w) {
                    distance[v] = distance[u] + w;
                    heap.update(Impl::nodeIndex(v));
                }
            });
    }
}

/**
 * Performs a Dijkstra traversal starting from the given source node.
 *
 * @param G The input graph.
 * @param u The source node.
 * @param lambda Takes a node and optionally its distance from the source node.
 * @throws std::runtime_error If the source node is not in the graph.
 */
template <typename Lambda, class GraphT>
    requires Impl::TraversalHandle<Lambda, typename GraphT::NodeT, typename GraphT::EdgeWeightT>
void DijkstraFrom(const GraphT &G, typename GraphT::NodeT u, Lambda lambda) {
    const auto sources = std::to_array({u});
    DijkstraFrom(G, sources.begin(), sources.end(), lambda);
}

} // namespace Traversal

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_DIJKSTRA_HPP_
