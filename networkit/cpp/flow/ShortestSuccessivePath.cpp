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
        [&](node u, node v, costs w, edgeid eid) { forwardEdges.emplace_back(u, v, w, eid); });

    // 3) Initial nodePotentials via Bellman–Ford on forward edges
    std::vector<double> nodePotentials(n, 0.0);
    for (index iter = 0; iter < n - 1; ++iter) {
        bool updated = false;
        for (auto &tuple : forwardEdges) {
            node u, v;
            edgeid e;
            double w;
            std::tie(u, v, e, w) = tuple;
            if (nodePotentials[u] + w < nodePotentials[v]) {
                nodePotentials[v] = nodePotentials[u] + w;
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
    std::vector<double> dist(n);
    std::vector<node> parentNode(n);
    std::vector<edgeid> parentEdge(n);
    std::vector<int> parentDir(n); // +1=forward, -1=backward

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
        std::fill(dist.begin(), dist.end(), infinity);
        dist[s] = 0;
        queue = decltype(queue)(); // clear
        queue.push({0, s});

        while (!queue.empty()) {
            auto [d, u] = queue.top();
            queue.pop();
            if (d > dist[u])
                continue;

            // — forward residual arcs u->v
            residualGraph.forEdgesOf(u, [&](node, node v, costs w, edgeid e) {
                double capRes = capAttr.get(e) - flowAttr.get(e);
                if (capRes <= epsilon)
                    return;
                double rawCost = w;
                double rc = rawCost + nodePotentials[u] - nodePotentials[v];
                if (dist[v] > dist[u] + rc) {
                    dist[v] = dist[u] + rc;
                    parentNode[v] = u;
                    parentEdge[v] = e;
                    parentDir[v] = +1;
                    queue.push({dist[v], v});
                }
            });

            // — backward residual arcs v->u
            residualGraph.forInEdgesOf(u, [&](node, node v, costs w, edgeid e) {
                double capRes = flowAttr.get(e);
                if (capRes <= epsilon)
                    return;
                double rawCost = -w;
                double rc = rawCost + nodePotentials[u] - nodePotentials[v];
                if (dist[v] > dist[u] + rc) {
                    dist[v] = dist[u] + rc;
                    parentNode[v] = u;
                    parentEdge[v] = e;
                    parentDir[v] = -1;
                    queue.push({dist[v], v});
                }
            });
        }

        // (c) update nodePotentials[u] += dist[u]
        for (node u = 0; u < n; ++u) {
            if (dist[u] < infinity) {
                nodePotentials[u] += dist[u];
            }
        }

        // (d) pick a demand node t reachable from s
        node t = none;
        for (node u = 0; u < n; ++u) {
            if (b[u] < -epsilon && dist[u] < infinity) {
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
                (parentDir[v] > 0) ? (capAttr.get(e) - flowAttr.get(e)) : flowAttr.get(e);
            f = std::min(f, capRes);
        }

        // (f) augment flow along the path
        for (node v = t; v != s; v = parentNode[v]) {
            edgeid e = parentEdge[v];
            double old = flowAttr.get(e);
            flowAttr.set(e, old + (parentDir[v] > 0 ? +f : -f));
        }

        // (g) update imbalances
        b[s] -= f;
        b[t] += f;
    }

    // At exit, flowAttr.get(eid) on each original edge is the computed min‐cost flow.
}

} // namespace NetworKit
