/*  ShortestSuccessivePath.hpp
 *
 *	Created on: 05.08.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/flow/ShortestSuccessivePath.hpp>

namespace NetworKit {

MinFlowShortestSuccessivePath::MinFlowShortestSuccessivePath(const Graph &G,
                                                             const std::string &capacityName,
                                                             const std::string &supplyName)
    : capacityAttributeName(capacityName), supplyAttributeName(supplyName), graph(&G) {

    if (!G.isDirected()) {
        throw std::runtime_error("MinFlowShortestSuccessivePath: Graph must be directed.");
    }

    if (!G.isWeighted()) {
        throw std::runtime_error("MinFlowShortestSuccessivePath: Graph must be weighted.");
    }

    if (!G.hasEdgeIds()) {
        throw std::runtime_error("MinFlowShortestSuccessivePath: Graph edges must be indexed.");
    }

    try {
        (void)G.edgeAttributes().find(capacityName);
    } catch (const std::runtime_error &e) {
        throw std::runtime_error("MinFlowShortestSuccessivePath: Provided edge attribute '"
                                 + capacityName + "' not found.");
    }

    try {
        (void)G.nodeAttributes().find(supplyName);
    } catch (const std::runtime_error &e) {
        throw std::runtime_error("MinFlowShortestSuccessivePath: Provided node attribute '"
                                 + supplyName + "' not found.");
    }
    residualGraph = *graph;
    auto flow = residualGraph.attachEdgeDoubleAttribute(FLOW);
    auto capacities = residualGraph.getEdgeDoubleAttribute(capacityName);
    auto supply = residualGraph.getNodeDoubleAttribute(supplyName);
    double totalSupply = 0.0;
    residualGraph.forNodes([&](node u) {
        totalSupply += supply.get(u);
        residualGraph.forEdgesOf(u, [&](node, node, cost, edgeid eid) {
            if (capacities.get(eid) < 0.0) {
                throw std::runtime_error(
                    "MinFlowShortestSuccessivePath: Capacities must be non-negative.");
            }
            flow.set(eid, 0.0);
        });
    });

    if (!Aux::NumericTools::equal(totalSupply, 0.0)) {
        throw std::runtime_error("MinFlowShortestSuccessivePath: Sum of node supplies and demands "
                                 "does not add up to zero.");
    }
}

const Graph::EdgeDoubleAttribute MinFlowShortestSuccessivePath::getFlow() const {
    if (!hasRun) {
        throw std::runtime_error(
            "MinFlowShortestSuccessivePath::getFlow: run() must be called first.");
    }
    return flows;
}

void MinFlowShortestSuccessivePath::run() {
    const count numberOfNodes = residualGraph.numberOfNodes();
    constexpr cost infiniteCosts = std::numeric_limits<cost>::infinity();
    constexpr double epsilon = 1e-12;

    // 1) Grab the three attributes:
    const auto capacities = residualGraph.getEdgeDoubleAttribute(capacityAttributeName);
    flows = residualGraph.getEdgeDoubleAttribute(FLOW);
    auto supply = residualGraph.getNodeDoubleAttribute(supplyAttributeName);

    // Apply Bellman-Ford to work to compute node potentials/distances
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
            "MinFlowShortestSuccessivePath: negative-cost cycle in residual graph");
    }

    // 5) Prepare Dijkstra data structures
    std::vector<cost> distances(numberOfNodes);
    std::vector<node> parentNode(numberOfNodes);
    std::vector<edgeid> parentEdge(numberOfNodes);
    std::vector<int> parentDirection(numberOfNodes); // +1=forward, -1=backward

    using costNodePair = std::pair<cost, node>;
    std::priority_queue<costNodePair, std::vector<costNodePair>, std::greater<costNodePair>> queue;

    // 6) Main successive‐shortest‐path loop
    while (true) {
        // (a) find a (new) supply node s with non-zero supply-value
        node s = none;
        for (node u = 0; u < numberOfNodes; ++u) {
            if (supply[u] > epsilon) {
                s = u;
                break;
            }
        }
        if (s == none)
            break;

        // (b) Dijkstra on residual network from s
        std::fill(distances.begin(), distances.end(), infiniteCosts);
        distances[s] = 0;
        queue = decltype(queue)(); // clear
        queue.emplace(0, s);
        while (!queue.empty()) {
            auto [distance, u] = queue.top();
            queue.pop();
            if (distance > distances[u])
                continue;

            // — forward residual arcs u->v
            residualGraph.forEdgesOf(u, [&](node, node v, cost cost, edgeid id) {
                const double residualCapacity = capacities.get(id) - flows.get(id);
                if (residualCapacity <= epsilon)
                    return;
                // The corrected costs are using nodePotentials to shift negative edges/costs to
                // positive ones
                const double correctedCost = cost + nodePotential[u] - nodePotential[v];
                if (distances[v] > distances[u] + correctedCost) {
                    distances[v] = distances[u] + correctedCost;
                    parentNode[v] = u;
                    parentEdge[v] = id;
                    parentDirection[v] = +1;
                    queue.emplace(distances[v], v);
                }
            });

            residualGraph.forInEdgesOf(u, [&](node, node v, cost cost, edgeid id) {
                double residualFlow = flows.get(id);
                if (residualFlow <= epsilon)
                    return;
                // The corrected costs are using nodePotentials to shift negative edges/costs to
                // positive ones
                const double correctedCost = -cost + nodePotential[u] - nodePotential[v];
                if (distances[v] > distances[u] + correctedCost) {
                    distances[v] = distances[u] + correctedCost;
                    parentNode[v] = u;
                    parentEdge[v] = id;
                    parentDirection[v] = -1;
                    queue.emplace(distances[v], v);
                }
            });
        }

        // (c) update nodePotentials[u] += dist[u]
        for (node u = 0; u < numberOfNodes; ++u) {
            if (distances[u] < infiniteCosts) {
                nodePotential[u] += distances[u];
            }
        }

        // (d) pick a demand node t reachable from s
        node t = none;
        for (node u = 0; u < numberOfNodes; ++u) {
            if (supply.get(u) < -epsilon && distances[u] < infiniteCosts) {
                t = u;
                break;
            }
        }

        if (t == none) {
            throw std::runtime_error(
                "MinFlowShortestSuccessivePath: unable to satisfy all supplies/demands");
        }

        // (e) compute bottleneck f = min(b[s], -b[t], min residual capacity on path)
        double f = std::min(supply.get(s), -supply.get(t));
        for (node v = t; v != s; v = parentNode[v]) {
            const edgeid e = parentEdge[v];
            const double capRes =
                (parentDirection[v] > 0) ? (capacities.get(e) - flows.get(e)) : flows.get(e);
            f = std::min(f, capRes);
        }

        // (f) augment flow along the path
        for (node v = t; v != s; v = parentNode[v]) {
            const edgeid e = parentEdge[v];
            const double old = flows.get(e);
            flows.set(e, old + (parentDirection[v] > 0 ? +f : -f));
        }

        // (g) update imbalances
        supply.set(s, supply.get(s) - f);
        supply.set(t, supply.get(t) + f);
    }
    totalCost = 0.0;

    residualGraph.forEdges(
        [&](node u, node v, cost c, edgeid eid) { totalCost += flows.get(eid) * c; });
    hasRun = true;
}

} // namespace NetworKit
