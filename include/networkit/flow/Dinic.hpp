/*  Dinic.hpp
 *
 *	Created on: 20.06.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_FLOW_DINIC_HPP_
#define NETWORKIT_FLOW_DINIC_HPP_
#include <vector>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
/**
 * @class Dinic
 * @brief Computes maximum flow in a directed, weighted graph using Dinic's algorithm.
 *
 * This class implements the blocking flow approach of Dinic's algorithm to
 * compute the maximum flow from a given source to a target node.
 */
class Dinic final : public Algorithm {

public:
    /**
     * @brief Construct a Dinic flow solver.
     *
     * @param G     Reference to the input directed, weighted graph.
     * @param src     Source node identifier.
     * @param dst     Target node identifier.
     * @throws std::runtime_error if the graph is not directed, not weighted, or src == dst.
     */
    Dinic(const Graph &G, node src, node dst);

    /**
     * @brief Execute the algorithm to compute maximum flow.
     *
     * Initializes the residual graph, then repeatedly performs BFS to build level
     * graphs and DFS to find blocking flows until no augmenting path remains.
     */
    void run() override;

    /**
     * @brief Retrieve the maximum flow value.
     *
     * @return The computed maximum flow from source to target.
     * @throws std::runtime_error if called before run().
     */
    edgeweight getMaxFlow() const;

private:
    void initializeResidualGraph();
    bool canReachTargetInLevelGraph();
    edgeweight computeBlockingPath();
    const Graph *graph;
    node source;
    node target;
    edgeweight maxFlow{};
    Graph residualGraph;
    std::vector<std::deque<node>> parents;
    static constexpr double RELATIVE_TOLERANCE = 1e-12;
    static constexpr double ABSOLUTE_TOLERANCE = 1e-15;
    edgeweight tolerance{};
};
} // namespace NetworKit
#endif // NETWORKIT_FLOW_DINIC_HPP_
