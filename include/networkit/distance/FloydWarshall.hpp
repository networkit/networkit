/*  FloydWarshall.hpp
 *
 *	Created on: 15.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_DISTANCE_FLOYD_WARSHALL_HPP_
#define NETWORKIT_DISTANCE_FLOYD_WARSHALL_HPP_

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @class FloydWarshall
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
template <class GraphType>
class FloydWarshall : public Algorithm {
    using NodeType = typename GraphType::node_type;
    using EdgeWeightType = typename GraphType::edge_weight_type;
    static constexpr NodeType null_node = std::numeric_limits<NodeType>::max();

public:
    /**
     * @brief Initializes the Floyd-Warshall algorithm for a given graph.
     *
     * The input graph must be weighted and may be either directed or undirected.
     *
     * @param G The graph on which shortest paths will be computed.
     */
    FloydWarshall(const GraphType &G) : graph(&G) {
        if (!G.isWeighted()) {
            throw std::runtime_error("The input graph is unweighted!");
        }
    }
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
    EdgeWeightType getDistance(NodeType source, NodeType target) const {
        assureFinished();
        return distances[source][target];
    }

    /**
     * @brief Checks whether a node is part of a negative cycle.
     *
     * A node is considered part of a negative cycle if its shortest path distance
     * to itself is negative.
     *
     * @param u The node to check.
     * @return `true` if the node is in a negative cycle, otherwise `false`.
     */
    bool isNodeInNegativeCycle(NodeType u) const {
        assureFinished();
        return nodesInNegativeCycle[u] == 1;
    }

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
    std::vector<NodeType> getNodesOnShortestPath(NodeType source, NodeType target) const;

private:
    const GraphType *graph;
    static constexpr EdgeWeightType infiniteDistance = std::numeric_limits<EdgeWeightType>::max();
    std::vector<std::vector<EdgeWeightType>> distances;
    std::vector<bool> nodesInNegativeCycle;
    std::vector<std::vector<NodeType>> pathMatrix;
    std::vector<std::vector<count>> hops;
    void tagNegativeCycles();
};

template <class GraphType>
void FloydWarshall<GraphType>::tagNegativeCycles() {
    graph->forNodes([&](NodeType w) {
        if (distances[w][w] >= 0.0)
            return;
        nodesInNegativeCycle[w] = 1;
        graph->forNodes([&](NodeType u) {
            if (distances[u][w] == infiniteDistance)
                return;
            graph->forNodes([&](NodeType v) {
                if (distances[w][v] != infiniteDistance) {
                    nodesInNegativeCycle[u] = 1;
                    nodesInNegativeCycle[v] = 1;
                    distances[u][v] = -std::numeric_limits<EdgeWeightType>::infinity();
                    pathMatrix[u][v] = null_node;
                }
            });
        });
    });
}

template <class GraphType>
void FloydWarshall<GraphType>::run() {
    const index numberOfNodes = graph->numberOfNodes();
    distances.resize(numberOfNodes, std::vector<EdgeWeightType>(numberOfNodes, infiniteDistance));
    nodesInNegativeCycle.resize(numberOfNodes);
    pathMatrix.resize(numberOfNodes, std::vector(numberOfNodes, null_node));
    hops.resize(numberOfNodes, std::vector(numberOfNodes, none));

    graph->forNodes([&](NodeType u) {
        distances[u][u] = 0.0;
        pathMatrix[u][u] = u;
        hops[u][u] = 0;
    });

    graph->forNodes([&](NodeType u) {
        for (const NodeType v : graph->neighborRange(u)) {
            distances[u][v] = graph->weight(u, v);
            pathMatrix[u][v] = v;
            hops[u][v] = 1;
        }
    });

    graph->forNodes([&](NodeType intermediate) {
        graph->parallelForNodes([&](NodeType source) {
            if (distances[source][intermediate] == infiniteDistance)
                return;
            graph->forNodes([&](NodeType target) {
                if (distances[intermediate][target] == infiniteDistance) {
                    return;
                }
                const EdgeWeightType candidateDistance =
                    distances[source][intermediate] + distances[intermediate][target];
                const count candidateHops = hops[source][intermediate] + hops[intermediate][target];
                if (candidateDistance < distances[source][target]) {
                    distances[source][target] = candidateDistance;
                    hops[source][target] = candidateHops;
                    pathMatrix[source][target] = pathMatrix[source][intermediate];
                }
                if (Aux::NumericTools::equal(candidateDistance, distances[source][target])
                    && candidateHops < hops[source][target]) {
                    hops[source][target] = candidateHops;
                    pathMatrix[source][target] = pathMatrix[source][intermediate];
                }
            });
        });
    });

    tagNegativeCycles();
    hasRun = true;
}

template <class GraphType>
std::vector<typename GraphType::node_type>
FloydWarshall<GraphType>::getNodesOnShortestPath(NodeType source, NodeType target) const {
    assureFinished();
    if (pathMatrix[source][target] == null_node) {
        return {};
    }
    std::vector<NodeType> path;
    NodeType currentNode = source;

    while (currentNode != target) {
        if (currentNode == null_node)
            return {};
        path.push_back(currentNode);
        currentNode = pathMatrix[currentNode][target];
    }
    path.push_back(target);
    return path;
}
} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_FLOYD_WARSHALL_HPP_
