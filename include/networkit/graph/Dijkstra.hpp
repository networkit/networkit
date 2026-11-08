
#ifndef NETWORKIT_GRAPH_DIJKSTRA_HPP_
#define NETWORKIT_GRAPH_DIJKSTRA_HPP_

#include <limits>
#include <vector>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/Graph.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

namespace NetworKit {

namespace Traversal {

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

/**
 * Iterate over nodes with Dijkstra starting from the nodes within the given range.
 *
 * @param G The input graph.
 * @param first The first element of the range.
 * @param last The end of the range.
 * @param handle Takes a node and its distance from the nodes in the range as input parameters.
 */
template <class InputIt, typename Handle, class GraphT>
void DijkstraFrom(const GraphT &G, InputIt first, InputIt last, Handle handle) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
    std::vector<EdgeWeightT> distance(G.upperNodeIdBound(),
                                      std::numeric_limits<EdgeWeightT>::max());
    const auto compareDistance = [&distance](NodeT u, NodeT v) noexcept -> bool {
        return distance[u] < distance[v];
    };

    tlx::d_ary_addressable_int_heap<NodeT, 2, decltype(compareDistance)> heap{compareDistance};

    std::for_each(first, last, [&heap, &distance](const NodeT u) {
        heap.push(u);
        distance[u] = 0;
    });

    while (!heap.empty()) {
        const auto u = heap.extract_top();
        callDijkstraHandle(handle, u, distance[u]);
        G.forNeighborsOf(
            u, [&](NodeT /* source */, const NodeT v, const EdgeWeightT w, edgeid /*id*/) {
                if (distance[v] > distance[u] + w) {
                    distance[v] = distance[u] + w;
                    heap.update(v);
                }
            });
    }
}

/**
 * Iterate over nodes with Dijkstra starting from the nodes within the given range.
 *
 * @param G The input graph.
 * @param u The source node.
 * @param lambda Takes a node and its distance from `u` as input parameters.
 */
template <typename Lambda, class GraphT>
void DijkstraFrom(const GraphT &G, typename GraphT::NodeT u, Lambda lambda) {
    using NodeT = typename GraphT::NodeT;
    std::vector<NodeT> vec({u});
    DijkstraFrom(G, vec.begin(), vec.end(), lambda);
}

} // namespace Traversal

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_DIJKSTRA_HPP_
