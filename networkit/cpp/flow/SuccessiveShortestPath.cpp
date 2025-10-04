/*  ShortestSuccessivePath.hpp
 *
 *	Created on: 05.08.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include "tlx/container/d_ary_heap.hpp"

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/flow/SuccessiveShortestPath.hpp>

namespace NetworKit {

SuccessiveShortestPathMinCostFlow::SuccessiveShortestPathMinCostFlow(const Graph &G,
                                                                     std::string_view capacityName,
                                                                     std::string_view supplyName)
    : graph(&G), capacityAttributeName(capacityName), supplyAttributeName(supplyName),
      totalCost(0) {

    if (!G.isDirected()) {
        throw std::runtime_error("SuccessiveShortestPathMinCostFlow: Graph must be directed.");
    }

    if (!G.isWeighted()) {
        throw std::runtime_error("SuccessiveShortestPathMinCostFlow: Graph must be weighted.");
    }

    if (!G.hasEdgeIds()) {
        throw std::runtime_error("SuccessiveShortestPathMinCostFlow: Graph edges must be indexed.");
    }

    try {
        (void)G.edgeAttributes().find(capacityAttributeName);
    } catch (const std::runtime_error &e) {
        throw std::runtime_error("SuccessiveShortestPathMinCostFlow: Provided edge attribute '"
                                 + capacityAttributeName + "' not found.");
    }

    try {
        (void)G.nodeAttributes().find(supplyAttributeName);
    } catch (const std::runtime_error &e) {
        throw std::runtime_error("SuccessiveShortestPathMinCostFlow: Provided node attribute '"
                                 + supplyAttributeName + "' not found.");
    }
    residualGraph = *graph;
    auto flow = residualGraph.attachEdgeDoubleAttribute(FLOW);
    auto capacities = residualGraph.getEdgeDoubleAttribute(capacityAttributeName);
    auto supply = residualGraph.getNodeDoubleAttribute(supplyAttributeName);
    double totalSupply = 0.0;
    residualGraph.forNodes([&](node u) {
        totalSupply += supply.get(u);
        residualGraph.forEdgesOf(u, [&](node, node, cost, edgeid eid) {
            if (capacities.get(eid) < 0.0) {
                throw std::runtime_error(
                    "SuccessiveShortestPathMinCostFlow: Capacities must be non-negative.");
            }
            flow.set(eid, 0.0);
        });
    });

    if (!Aux::NumericTools::equal(totalSupply, 0.0)) {
        throw std::runtime_error(
            "SuccessiveShortestPathMinCostFlow: Sum of node supplies and demands "
            "does not add up to zero.");
    }
}

std::vector<SuccessiveShortestPathMinCostFlow::cost>
SuccessiveShortestPathMinCostFlow::computeNodePotentials(count numberOfNodes) const {
    // Apply Bellman-Ford to compute node potentials/distances (dealing with negative weights/costs)
    std::vector<cost> nodePotential(numberOfNodes, 0.0);
    for (count i = 1; i < numberOfNodes; ++i) {
        bool updated = false;
        residualGraph.forEdges([&](node u, node v, cost cost, edgeid) {
            if (nodePotential[u] + cost < nodePotential[v]) {
                updated = true;
                nodePotential[v] = nodePotential[u] + cost;
            }
        });
        if (!updated)
            break;
    }

    // negative cycle detection
    // Note: We can not throw in parallel loop since it will terminate the program
    bool negativeCycleDetected = false;
    residualGraph.parallelForEdges([&](node u, node v, cost cost, edgeid) {
        if (nodePotential[u] + cost < nodePotential[v]) {
            negativeCycleDetected = true;
        }
    });

    if (negativeCycleDetected) {
        throw std::runtime_error(
            "SuccessiveShortestPathMinCostFlow: negative-cost cycle in residual graph");
    }
    return nodePotential;
}

void SuccessiveShortestPathMinCostFlow::dijkstraOnResidualGraph(
    node start, std::span<const cost> nodePotential, std::span<cost> distances,
    std::span<node> parentNode, std::span<edgeid> parentEdge, std::span<int> parentDirection,
    const Graph::EdgeDoubleAttribute &capacities) const {

    const auto comparator = [&](const NodeWithCost &nwc1, const NodeWithCost &nwc2) -> bool {
        if (!Aux::NumericTools::equal(nwc1.distance, nwc2.distance))
            return nwc1.distance > nwc2.distance;
        return nwc1.u > nwc2.u;
    };

    tlx::d_ary_heap<NodeWithCost, 2, decltype(comparator)> minHeap(comparator);

    std::ranges::fill(distances, infiniteCosts);
    distances[start] = 0;
    minHeap.clear();
    minHeap.push({0, start});
    do {
        const cost distance = minHeap.top().distance;
        const node u = minHeap.top().u;
        minHeap.pop();
        if (distance > distances[u])
            continue;

        // — forward residual arcs u->v
        residualGraph.forEdgesOf(u, [&](node, node v, cost cost, edgeid eid) {
            const double residualCapacity = capacities.get(eid) - flows.get(eid);
            if (residualCapacity <= epsilon)
                return;
            // The corrected costs are using nodePotentials to shift negative edges/costs to
            // positive ones
            const double correctedCost = cost + nodePotential[u] - nodePotential[v];
            if (distances[v] > distances[u] + correctedCost) {
                distances[v] = distances[u] + correctedCost;
                parentNode[v] = u;
                parentEdge[v] = eid;
                parentDirection[v] = +1;
                minHeap.push({distances[v], v});
            }
        });

        // — backward residual arcs v->u
        residualGraph.forInEdgesOf(u, [&](node, node v, cost cost, edgeid eid) {
            double residualFlow = flows.get(eid);
            if (residualFlow <= epsilon)
                return;
            // The corrected costs are using nodePotentials to shift negative edges/costs to
            // positive ones
            const double correctedCost = -cost + nodePotential[u] - nodePotential[v];
            if (distances[v] > distances[u] + correctedCost) {
                distances[v] = distances[u] + correctedCost;
                parentNode[v] = u;
                parentEdge[v] = eid;
                parentDirection[v] = -1;
                minHeap.push({distances[v], v});
            }
        });
    } while (!minHeap.empty());
}

void SuccessiveShortestPathMinCostFlow::run() {
    const count numberOfNodes = residualGraph.numberOfNodes();

    const auto capacities = residualGraph.getEdgeDoubleAttribute(capacityAttributeName);
    flows = residualGraph.getEdgeDoubleAttribute(FLOW);
    auto supply = residualGraph.getNodeDoubleAttribute(supplyAttributeName);

    std::vector<cost> nodePotential = computeNodePotentials(numberOfNodes);

    // Prepare Dijkstra data structures
    std::vector<cost> distances(numberOfNodes);
    std::vector<node> parentNode(numberOfNodes);
    std::vector<edgeid> parentEdge(numberOfNodes);
    std::vector<int> parentDirection(numberOfNodes);

    // Main successive‐shortest‐path loop
    do {
        // (a) find a (new) supply node start with non-zero supply-value
        node start = none;
        for (const node u : residualGraph.nodeRange()) {
            if (supply[u] > epsilon) {
                start = u;
                break;
            }
        }
        if (start == none)
            break;

        // (b) Dijkstra on residual network from start node
        dijkstraOnResidualGraph(start, nodePotential, distances, parentNode, parentEdge,
                                parentDirection, capacities);

        // (c) update nodePotentials
        residualGraph.parallelForNodes([&](node u) {
            if (distances[u] < infiniteCosts) {
                nodePotential[u] += distances[u];
            }
        });

        // (d) pick a demand node target reachable from start
        node target = none;
        for (const node u : residualGraph.nodeRange()) {
            if (supply.get(u) < -epsilon && distances[u] < infiniteCosts) {
                target = u;
                break;
            }
        }

        if (target == none) {
            throw std::runtime_error(
                "SuccessiveShortestPathMinCostFlow: unable to satisfy all supplies/demands");
        }

        // (e) compute bottleneck-flow = min(b[s], -b[t], min residual capacity on path)
        double bottleneckFlow = std::min(supply.get(start), -supply.get(target));
        for (node v = target; v != start; v = parentNode[v]) {
            const edgeid e = parentEdge[v];
            const double residualCapacity =
                (parentDirection[v] > 0) ? (capacities.get(e) - flows.get(e)) : flows.get(e);
            bottleneckFlow = std::min(bottleneckFlow, residualCapacity);
        }

        // (f) augment flow along the path
        for (node v = target; v != start; v = parentNode[v]) {
            const edgeid e = parentEdge[v];
            const double old = flows.get(e);
            flows.set(e, old + (parentDirection[v] > 0 ? +bottleneckFlow : -bottleneckFlow));
        }

        // (g) update imbalances
        supply.set(start, supply.get(start) - bottleneckFlow);
        supply.set(target, supply.get(target) + bottleneckFlow);
    } while (true);

    totalCost = residualGraph.parallelSumForEdges(
        [&](node /*u*/, node /*v*/, cost c, edgeid eid) -> cost { return flows.get(eid) * c; });
    hasRun = true;
}

} // namespace NetworKit
