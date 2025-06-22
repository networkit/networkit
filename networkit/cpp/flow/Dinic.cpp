/*  Dinic.cpp
 *
 *	Created on: 20.06.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <deque>
#include <queue>
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

    residualGraph = Graph(graph->upperNodeIdBound(), true, true);
}

void Dinic::buildResidual() {
    graph->forEdges([&](node u, node v, edgeweight w) {
        residualGraph.addEdge(u, v, w);
        residualGraph.addEdge(v, u, 0.0);
    });
}

bool Dinic::get_parents_bfs() {
    std::fill(levels.begin(), levels.end(), -1);
    for (auto &deque : parents) {
        deque.clear();
    }
    std::queue<node> queue;
    levels[source] = 0;
    queue.push(source);
    while (!queue.empty()) {
        const node u = queue.front();
        queue.pop();
        for (const node v : residualGraph.neighborRange(u)) {
            if (residualGraph.weight(u, v) > 0.0) {
                if (levels[v] == -1) {
                    levels[v] = levels[u] + 1;
                    parents[v].push_back(u);
                    queue.push(v);
                } else if (levels[v] == levels[u] + 1) {
                    parents[v].push_back(u);
                }
            }
        }
    }
    return levels[target] >= 0;
}

edgeweight Dinic::blocking_path() {
    edgeweight totalFlow = 0.0;
    std::vector<node> path;
    path.reserve(residualGraph.numberOfNodes());
    path.push_back(target);
    node u = target;
    while (true) {
        node v;
        if (!parents[u].empty()) {
            v = parents[u].front();
            path.push_back(v);
        } else {
            path.pop_back();
            if (path.empty())
                break;
            v = path.back();
            parents[v].pop_front();
        }

        // path build from target to source
        if (v == source) {
            edgeweight minimalFlowOnPath = std::numeric_limits<edgeweight>::max();
            for (int i{}; i + 1 < path.size(); ++i) {
                const node parent = path[i];
                const node child = path[i + 1];
                minimalFlowOnPath =
                    std::min(minimalFlowOnPath, residualGraph.weight(parent, child));
            }
            for (int i{}; i + 1 < path.size(); ++i) {
                const node parent = path[i];
                const node child = path[i + 1];
                const edgeweight currentCapacity = residualGraph.weight(parent, child);
                residualGraph.setWeight(parent, child, currentCapacity - minimalFlowOnPath);
                if (residualGraph.hasEdge(child, parent)) {
                    const edgeweight reverseCapacity = residualGraph.weight(child, parent);
                    residualGraph.setWeight(child, parent, reverseCapacity + minimalFlowOnPath);
                } else {
                    residualGraph.addEdge(child, parent, minimalFlowOnPath);
                }
                if (residualGraph.weight(parent, child) == 0) {
                    parents[child].pop_front();
                }
            }
            totalFlow += minimalFlowOnPath;
            v = path.back();
        }
        u = v;
    }
    return totalFlow;
}

void Dinic::run() {
    buildResidual();
    maxFlow = 0.0;
    while (get_parents_bfs()) {
        maxFlow += blocking_path();
    }
    hasRun = true;
}

} // namespace NetworKit
