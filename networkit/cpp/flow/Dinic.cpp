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

namespace NetworKit {

Dinic::Dinic(const Graph &G, node s, node t) : graph(&G), source(s), target(t) {

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
    residualGraph = Graph(graph->upperNodeIdBound(), true, true);
}

void Dinic::initializeResidualGraph() {
    graph->forEdges([&](node u, node v, edgeweight w) {
        // add original capacity
        residualGraph.addEdge(u, v, w);
        // add reverse edge to track flow, initially 0
        residualGraph.addEdge(v, u, 0.0);
    });
}

bool Dinic::determineValidParents() {
    std::vector<int> level(residualGraph.numberOfNodes(), -1);
    for (auto &parentList : parents) {
        parentList.clear();
    }
    std::queue<node> queue;
    level[source] = 0;
    queue.push(source);
    while (!queue.empty()) {
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
    }
    return level[target] > -1;
}

edgeweight Dinic::computeBlockingPath() {
    edgeweight totalFlow = 0.0;
    std::vector<node> path;
    path.push_back(target);
    node u = target;
    while (true) {
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
            for (size_t i{}; i + 1 < path.size(); ++i) {
                const node parent = path[i + 1];
                const node child = path[i];
                bottleNeckOnPath = std::min(bottleNeckOnPath, residualGraph.weight(parent, child));
            }
            // update the capacities and flows in the other edges
            for (size_t i{}; i + 1 < path.size(); ++i) {
                const node parent = path[i + 1];
                const node child = path[i];
                const edgeweight currentCapacity = residualGraph.weight(parent, child);
                residualGraph.setWeight(parent, child, currentCapacity - bottleNeckOnPath);
                if (residualGraph.hasEdge(child, parent)) {
                    const edgeweight reverseCapacity = residualGraph.weight(child, parent);
                    residualGraph.setWeight(child, parent, reverseCapacity + bottleNeckOnPath);
                } else {
                    residualGraph.addEdge(child, parent, bottleNeckOnPath);
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
    }

    return totalFlow;
}

void Dinic::run() {
    initializeResidualGraph();
    maxFlow = 0.0;
    while (determineValidParents()) {
        if (const double flow = computeBlockingPath(); !Aux::NumericTools::equal(flow, 0.0)) {
            maxFlow += flow;
        } else
            break;
    }
    hasRun = true;
}

edgeweight Dinic::getMaxFlow() const {
    assureFinished();
    return maxFlow;
}

} // namespace NetworKit
