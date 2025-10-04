/*  ShortestSuccessivePath.hpp
 *
 *	Created on: 05.08.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_FLOW_SUCCESSIVE_SHORTEST_PATH_HPP_
#define NETWORKIT_FLOW_SUCCESSIVE_SHORTEST_PATH_HPP_

#include <span>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
class SuccessiveShortestPathMinCostFlow : public Algorithm {
    using cost = edgeweight;

public:
    /**
     * @brief Constructs a min-cost flow solver using the successive shortest path algorithm.
     *
     * Initializes the solver on a given directed, weighted graph with edge capacities
     * and node supply/demand attributes.
     *
     * @param G             Underlying graph (must be directed, weighted, and have indexed edges).
     * @param capacityName  Name of the edge attribute holding non-negative capacity values.
     * @param supplyName    Name of the node attribute holding supply/demand values.
     *
     * @throws std::runtime_error if:
     *         - the graph is not directed,
     *         - the graph is not weighted,
     *         - the graph does not have indexed edges,
     *         - the specified edge attribute (capacities) does not exist,
     *         - the specified node attribute (supplies) does not exist,
     *         - any edge capacity is negative,
     *         - or the sum of all node supplies/demands is not zero.
     */
    SuccessiveShortestPathMinCostFlow(const Graph &G, std::string_view capacityName,
                                      std::string_view supplyName);

    /**
     * \brief Compute a minimum-cost feasible flow using the successive shortest path algorithm.
     *
     * \details
     * This method implements the classic Successive Shortest Path (SSP) algorithm:
     *  - Initializes node potentials using Bellman–Ford to reweight edges
     *    (ensuring non-negative reduced costs).
     *  - Iteratively selects a supply node with positive imbalance.
     *  - Runs Dijkstra on the residual network (with reduced costs) to find
     *    the cheapest augmenting path to a reachable demand node.
     *  - Augments as much flow as possible along that path (bounded by supply,
     *    demand, and residual capacities).
     *  - Updates node potentials and node imbalances.
     *
     * The process continues until all supplies are pushed to demands.
     * At the end, the edge attribute \c FLOW contains the computed flow on
     * each edge, and \c getTotalCost() can be used to retrieve the total cost.
     *
     * \throws std::runtime_error if:
     *  - a negative-cost cycle is detected in the residual graph,
     *  - or if there exists remaining supply/demand that cannot be routed to any demand/supply.
     */
    void run() override;

    /**
     * @brief Returns the total cost of the computed minimum-cost flow.
     *
     * @return Total flow cost as a double.
     *
     * @throws std::runtime_error if called before run() has finished.
     */
    double getTotalCost() const {
        assureFinished();
        return totalCost;
    }

    /**
     * @brief Returns the flow values assigned to each edge.
     *
     * The returned attribute maps every edge to the computed flow value.
     * Flow values are always non-negative and do not exceed the specified capacity.
     *
     * @return EdgeDoubleAttribute object with flow values per edge.
     *
     * @throws std::runtime_error if called before run() has finished.
     */
    Graph::EdgeDoubleAttribute getFlow() const {
        assureFinished();
        return flows;
    }

private:
    struct NodeWithCost {
        cost distance;
        node u;
    };
    std::vector<cost> computeNodePotentials(count numberOfNodes) const;
    void dijkstraOnResidualGraph(node start, std::span<const cost> nodePotential,
                                 std::span<cost> distances, std::span<node> parentNode,
                                 std::span<edgeid> parentEdge, std::span<int> parentDirection,
                                 const Graph::EdgeDoubleAttribute &capacities) const;
    const Graph *graph;
    const std::string capacityAttributeName;
    const std::string supplyAttributeName;
    Graph residualGraph;
    static constexpr const char *FLOW = "flow";
    Graph::EdgeDoubleAttribute flows;
    cost totalCost;
    static constexpr cost infiniteCosts = std::numeric_limits<cost>::infinity();
    static constexpr double epsilon = 1e-12;
};

} // namespace NetworKit
#endif // NETWORKIT_FLOW_SUCCESSIVE_SHORTEST_PATH_HPP_
