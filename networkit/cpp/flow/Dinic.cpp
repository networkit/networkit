/*  Dinic.cpp
 *
 *	Created on: 20.06.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <deque>
#include <queue>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/flow/Dinic.hpp>
#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

Dinic::Dinic(const Graph &G, node src, node dst) : graph(&G), source(src), target(dst) {

    if (!graph->isDirected()) {
        throw std::runtime_error("Dinic algorithm requires directed graph!");
    }
    if (!graph->isWeighted()) {
        throw std::runtime_error("Dinic algorithm requires weighted graph!");
    }
    if (source == target) {
        throw std::runtime_error(
            "Dinic algorithm requires `source` and `target` node to be different!");
    }
    parents.resize(graph->numberOfNodes());
}

void Dinic::initializeResidualGraph() {
    // Start from the original graph.
    residualGraph = *graph;

    // Add missing reverse arcs with 0 capacity (but don't overwrite real antiparallel edges).
    graph->forEdges([&](node u, node v, edgeweight w) {
        if (w < 0.0) {
            throw std::runtime_error("Dinic requires non-negative capacities!");
        }
        if (!residualGraph.hasEdge(v, u)) {
            residualGraph.addEdge(v, u, 0.0);
        }
    });
    // Rebuild edge indices after structural changes.
    residualGraph.indexEdges();
}

bool Dinic::canReachTargetInLevelGraph() {
    std::vector<int> level(residualGraph.numberOfNodes(), -1);
    for (auto &parentList : parents) {
        parentList.clear();
    }
    std::queue<node> queue;
    level[source] = 0;
    queue.push(source);
    do {
        const node parent = queue.front();
        queue.pop();
        for (const node child : residualGraph.neighborRange(parent)) {
            // We only consider connections with positive remaining capacity
            if (residualGraph.weight(parent, child) > 0.0) {
                if (level[child] == -1) {
                    level[child] = level[parent] + 1;
                    parents[child].push_back(parent);
                    queue.push(child);
                } else if (level[child] == level[parent] + 1) {
                    parents[child].push_back(parent);
                }
            }
        }
    } while (!queue.empty());

    return level[target] > -1;
}

edgeweight Dinic::computeBlockingPath() {
    edgeweight totalFlow = 0.0;
    std::vector<node> path;
    path.push_back(target);
    node u = target;
    auto flow = residualGraph.getEdgeDoubleAttribute(FLOW);
    do {
        node v = none;
        // build path from target to source
        if (!parents[u].empty()) {
            v = parents[u].front();
            path.push_back(v);
        } else {
            path.pop_back();
            if (path.empty())
                break;
            v = path.back();
        }
        // path has been build from target to source, so the parent is on i+1 position of the ith
        // child
        if (v == source) {
            edgeweight bottleNeckOnPath = std::numeric_limits<edgeweight>::max();
            // determine minimal flow on path
            for (size_t i = 0; i + 1 < path.size(); ++i) {
                const node parent = path[i + 1];
                const node child = path[i];
                bottleNeckOnPath = std::min(bottleNeckOnPath, residualGraph.weight(parent, child));
            }
            // update the capacities and flows in the other edges
            for (size_t i = 0; i + 1 < path.size(); ++i) {
                const node parent = path[i + 1];
                const node child = path[i];
                const index edgeID = residualGraph.edgeId(parent, child);
                residualGraph.setWeight(parent, child,
                                        residualGraph.weight(parent, child) - bottleNeckOnPath);
                if (flow.get(edgeID) > 0.0) {
                    flow.set(edgeID, flow.get(edgeID) + bottleNeckOnPath);
                } else {
                    flow.set(edgeID, bottleNeckOnPath);
                }
                if (residualGraph.weight(parent, child) == 0 && !parents[child].empty()) {
                    parents[child].pop_front();
                }
            }
            totalFlow += bottleNeckOnPath;
            path.clear();
            path.push_back(target);
        }
        u = v;
    } while (true);

    return totalFlow;
}

void Dinic::run() {
    initializeResidualGraph();
    maxFlow = 0.0;
    while (canReachTargetInLevelGraph()) {
        const double flow = computeBlockingPath();
        if (Aux::NumericTools::equal(flow, 0.0))
            break;
        maxFlow += flow;
    }
    hasRun = true;
}

edgeweight Dinic::getMaxFlow() const {
    assureFinished();
    return maxFlow;
}

} // namespace NetworKit
