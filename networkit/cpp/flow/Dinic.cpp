/*  Dinic.cpp
 *
 *	Created on: 20.06.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <deque>
#include <queue>
#include <stack>
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
        throw std::runtime_error("Dinic algorithm requires `source` and `target` node to be different!");
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
    for (auto & deque : parents) {
        deque.clear();
    }
    std::queue<node> queue;
    levels[source] = 0;
    queue.push(source);
    while (!queue.empty()) {
        const node u = queue.front();
        queue.pop();
        for (const node v: residualGraph.neighborRange(u)) {
            if (residualGraph.weight(u, v) > 0.0) {
                if(levels[v] == -1) {
                    levels[v] = levels[u] + 1;
                    parents[v].push_back(u);
                    queue.push(v);
                }
                else if (levels[v] == levels[u] + 1 ) {
                    parents[v].push_back(u);
                }
            }
        }
    }
    return levels[target] >= 0;
}

edgeweight Dinic::dfs(node u, edgeweight flow) {
}

void Dinic::run() {
    buildResidual();
    edgeweight flow = 0.0;
    level.resize(residualGraph.upperNodeIdBound());
    ptr.resize(residualGraph.upperNodeIdBound());
    // Main loop: while there is an augmenting path
    while (bfs()) {
        std::ranges::fill(ptr, 0);
        // Send blocking flow
        edgeweight pushed;
        do {
            pushed = dfs(source, std::numeric_limits<edgeweight>::max());
            maxFlow += pushed;
        } while (pushed > 0.0);
    }
    maxFlow = flow;
    hasRun = true;
}
} // namespace NetworKit
