/*  ShortestSuccessivePath.hpp
 *
 *	Created on: 05.08.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <networkit/flow/ShortestSuccessivePath.hpp>
namespace NetworKit {

MinFlowShortestSuccessivePath::MinFlowShortestSuccessivePath(const Graph &G,
                                                             const std::string &capacityName,
                                                             const std::string &supplyName)
    : capacityAttributeName(capacityName), supplyAttributeName(supplyName), graph(&G) {

    if (!G.isDirected()) {
        throw std::runtime_error("MinFlowShortestSuccessivePath: Graph must be directed");
    }

    if (!G.isWeighted()) {
        throw std::runtime_error("MinFlowShortestSuccessivePath: Graph must be weighted.");
    }

    if (!G.hasEdgeIds()) {
        throw std::runtime_error("MinFlowShortestSuccessivePath: Graph edges must be indexed");
    }

    try {
        // this internally calls AttributeMap::find and will throw if the attribute does not exist
        (void)G.edgeAttributes().find(capacityName);
    } catch (const std::runtime_error &e) {
        // catch the generic “No such attribute” and give our own
        throw std::runtime_error("MinFlowShortestSuccessivePath: Provided edge attribute '"
                                 + capacityName + "' not found");
    }

    try {
        // this internally calls AttributeMap::find and will throw if the attribute does not exist
        (void)G.nodeAttributes().find(supplyName);
    } catch (const std::runtime_error &e) {
        // catch the generic “No such attribute” and give our own
        throw std::runtime_error("MinFlowShortestSuccessivePath: Provided node attribute '"
                                 + supplyName + "' not found");
    }
    residualGraph = *graph;
    auto flow = residualGraph.attachEdgeDoubleAttribute(FLOW);
    auto capacities = residualGraph.getEdgeDoubleAttribute(capacityName);
    residualGraph.forNodes([&](node u) {
        residualGraph.forEdgesOf(u, [&](node, node, costs, edgeid eid) {
            if (capacities.get(eid) < 0.0) {
                throw std::runtime_error(
                    "MinFlowShortestSuccessivePath: Capacities must be non-negative");
            }
            flow.set(eid, 0.0);
        });
    });
}

void MinFlowShortestSuccessivePath::run() {
    const index n = residualGraph.numberOfNodes();
    constexpr double infinity = std::numeric_limits<double>::infinity();
    constexpr double epsilon = 1e-12;

    // 1) Grab the three attributes:
    auto capAttr = residualGraph.getEdgeDoubleAttribute(capacityAttributeName);
    auto flowAttr = residualGraph.getEdgeDoubleAttribute(FLOW);
    auto supplyAttr = residualGraph.getNodeDoubleAttribute(supplyAttributeName);

    // 2) Gather all forward edges (u->v, eid, cost=w) for Bellman–Ford
    std::vector<std::tuple<node, node, costs, edgeid>> forwardEdges;
    residualGraph.forEdges(
        [&](node u, node v, costs cost, edgeid eid) { forwardEdges.emplace_back(u, v, cost, eid); });

    // 3) Initial nodePotentials via Bellman–Ford on forward edges
    std::vector<double> nodePotentials(n, 0.0);
    for (index iter = 0; iter < n - 1; ++iter) {
        bool updated = false;
        for (auto &tuple : forwardEdges) {
            node u, v;
            edgeid e;
            costs cost;
            std::tie(u, v, cost, e) = tuple;
            if (nodePotentials[u] + cost < nodePotentials[v]) {
                nodePotentials[v] = nodePotentials[u] + cost;
                updated = true;
            }
        }
        if (!updated)
            break;
    }

    // 4) Read node imbalances b[u] = supply[u]
    std::vector<double> b(n);
    residualGraph.forNodes([&](node u) { b[u] = supplyAttr.get(u); });

    // 5) Prepare Dijkstra data structures
    std::vector<double> distances(n);
    std::vector<node> parentNode(n);
    std::vector<edgeid> parentEdge(n);
    std::vector<int> parentDirection(n); // +1=forward, -1=backward

    using costNodePair = std::pair<costs, node>;
    std::priority_queue<costNodePair, std::vector<costNodePair>, std::greater<costNodePair>> queue;

    // 6) Main successive‐shortest‐path loop
    while (true) {
        // (a) find a supply node s
        node s = none;
        for (node u = 0; u < n; ++u) {
            if (b[u] > epsilon) {
                s = u;
                break;
            }
        }
        if (s == none)
            break; // done

        // (b) Dijkstra on residual network from s
        std::fill(distances.begin(), distances.end(), infinity);
        distances[s] = 0;
        queue = decltype(queue)(); // clear
        queue.push({0, s});
        while (!queue.empty()) {
            auto [distance, u] = queue.top();
            queue.pop();
            if (distance > distances[u])
                continue;


            // — forward residual arcs u->v
            residualGraph.forEdgesOf(u, [&](node, node v, costs cost, edgeid id) {
                const double residualCapacity = capAttr.get(id) - flowAttr.get(id);
                if (residualCapacity <= epsilon)
                    return;
                const double rawCost = cost + nodePotentials[u] - nodePotentials[v];
                if (distances[v] > distances[u] + rawCost) {
                    distances[v] = distances[u] + rawCost;
                    parentNode[v] = u;
                    parentEdge[v] = id;
                    parentDirection[v] = +1;
                    queue.emplace(distances[v], v);
                }
            });

            // — backward residual arcs v->u
            residualGraph.forInEdgesOf(u, [&](node, node v, costs cost, edgeid id) {
                double capRes = flowAttr.get(id);
                if (capRes <= epsilon)
                    return;
                double rawCost = -cost;
                double rc = rawCost + nodePotentials[u] - nodePotentials[v];
                if (distances[v] > distances[u] + rc) {
                    distances[v] = distances[u] + rc;
                    parentNode[v] = u;
                    parentEdge[v] = id;
                    parentDirection[v] = -1;
                    queue.push({distances[v], v});
                }
            });
        }

        // (c) update nodePotentials[u] += dist[u]
        for (node u = 0; u < n; ++u) {
            if (distances[u] < infinity) {
                nodePotentials[u] += distances[u];
            }
        }

        // (d) pick a demand node t reachable from s
        node t = none;
        for (node u = 0; u < n; ++u) {
            if (b[u] < -epsilon && distances[u] < infinity) {
                t = u;
                break;
            }
        }
        if (t == none) {
            throw std::runtime_error(
                "MinFlowShortestSuccessivePath: unable to satisfy all supplies");
        }

        // (e) compute bottleneck f = min(b[s], -b[t], min residual capacity on path)
        double f = std::min(b[s], -b[t]);
        for (node v = t; v != s; v = parentNode[v]) {
            edgeid e = parentEdge[v];
            double capRes =
                (parentDirection[v] > 0) ? (capAttr.get(e) - flowAttr.get(e)) : flowAttr.get(e);
            f = std::min(f, capRes);
        }

        // (f) augment flow along the path
        for (node v = t; v != s; v = parentNode[v]) {
            edgeid e = parentEdge[v];
            double old = flowAttr.get(e);
            flowAttr.set(e, old + (parentDirection[v] > 0 ? +f : -f));
        }

        // (g) update imbalances
        b[s] -= f;
        b[t] += f;
    }
    totalCost = 0.0;

    residualGraph.forEdges([&](node u, node v, costs c, edgeid eid) {
        double f = flowAttr.get(eid);
        totalCost += f * c;
    });
    hasRun = true;
}

} // namespace NetworKit
