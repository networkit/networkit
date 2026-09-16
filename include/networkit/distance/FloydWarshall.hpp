/*  FloydWarshall.hpp
 *
 *	Created on: 15.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_DISTANCE_FLOYD_WARSHALL_HPP_
#define NETWORKIT_DISTANCE_FLOYD_WARSHALL_HPP_

#include <limits>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @class GenericFloydWarshall
 * @brief Computes all-pairs shortest paths using the Floyd-Warshall algorithm.
 *
 * This algorithm finds the shortest paths between all node pairs in a weighted graph,
 * supporting both directed and undirected graphs. It correctly handles negative edge
 * weights and detects negative cycles. If multiple shortest paths exist, it returns
 * one with the fewest nodes.
 *
 * The algorithm has a time complexity of O(n³), making it suitable for small to
 * medium-sized graphs.
 */
template <typename GraphT>
class GenericFloydWarshall : public Algorithm {
public:
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
    using DistanceT = edgeweight;

    /**
     * @brief Initializes the Floyd-Warshall algorithm for a given graph.
     *
     * The input graph must be weighted and may be either directed or undirected.
     *
     * @param G The graph on which shortest paths will be computed.
     */
    GenericFloydWarshall(const GraphT &G);
    /**
     * @brief Runs the Floyd-Warshall algorithm.
     *
     * Computes shortest path distances and reconstructs paths between all node pairs.
     * Also identifies nodes involved in negative cycles.
     */
    void run() override;

    /**
     * @brief Returns the shortest distance between two nodes.
     *
     * If no path exists, returns `std::numeric_limits<edgeweight>::max()`.
     *
     * @param source The starting node.
     * @param target The destination node.
     * @return The shortest path distance from `source` to `target`.
     */
    DistanceT getDistance(NodeT source, NodeT target) const;
    /**
     * @brief Checks whether a node is part of a negative cycle.
     *
     * A node is considered part of a negative cycle if its shortest path distance
     * to itself is negative.
     *
     * @param u The node to check.
     * @return `true` if the node is in a negative cycle, otherwise `false`.
     */
    bool isNodeInNegativeCycle(NodeT u) const;

    /**
     * @brief Retrieves the shortest path between two nodes.
     *
     * Returns a sequence of nodes representing the shortest path from `source` to
     * `target`. If no path exists, returns an empty vector.
     *
     * If multiple shortest paths exist with the same total distance, the function
     * returns the one with the fewest nodes.
     *
     * @param source The starting node.
     * @param target The destination node.
     * @return A vector of nodes forming the shortest path.
     */
    std::vector<NodeT> getNodesOnShortestPath(NodeT source, NodeT target) const;

    /**
     * @brief Returns the full all-pairs distance matrix.
     *
     * The returned object is a reference to internal storage. Do not modify it
     * directly; call run() again if the graph changes.
     *
     * @return A matrix where entry [u][v] is the shortest distance from u to v.
     */
    const std::vector<std::vector<DistanceT>> &getDistances() const;

private:
    const GraphT *graph;
    static constexpr DistanceT infiniteDistance = std::numeric_limits<DistanceT>::max();
    std::vector<std::vector<DistanceT>> distances;
    std::vector<bool> nodesInNegativeCycle;
    std::vector<std::vector<NodeT>> pathMatrix;
    std::vector<std::vector<count>> hops;
    void tagNegativeCycles();
};

using FloydWarshall = GenericFloydWarshall<Graph>;

} // namespace NetworKit

#include <networkit/distance/FloydWarshallImpl.hpp>

#endif // NETWORKIT_DISTANCE_FLOYD_WARSHALL_HPP_
