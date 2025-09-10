/*  ShortestSuccessivePath.hpp
 *
 *	Created on: 05.08.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_FLOW_SHORTEST_SUCCESSIVE_PATH_HPP_
#define NETWORKIT_FLOW_SHORTEST_SUCCESSIVE_PATH_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
class MinFlowShortestSuccessivePath : public Algorithm {
    using cost = edgeweight;

public:
    /**
     * @param G             underlying graph (must be directed, weighted, and have indexed edges)
     * @param capacityName  name of the edge attribute holding non-negative capacity values
     * @param supplyName    name of the node attribute holding supply/demand values
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
    MinFlowShortestSuccessivePath(const Graph &G, const std::string &capacityName,
                                  const std::string &supplyName);

    void run() override;
    double getTotalCost() const {
        assureFinished();
        return totalCost;
    }
    const Graph::EdgeDoubleAttribute getFlow() const;

private:
    const Graph *graph;
    std::string capacityAttributeName;
    std::string supplyAttributeName;
    Graph residualGraph;
    static constexpr const char *FLOW = "flow";
    Graph::EdgeDoubleAttribute flows;
    cost totalCost{};
};

} // namespace NetworKit
#endif // NETWORKIT_FLOW_SHORTEST_SUCCESSIVE_PATH_HPP_
