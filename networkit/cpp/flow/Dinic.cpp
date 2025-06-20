/*  Dinic.cpp
 *
 *	Created on: 20.06.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <queue>
#include <stack>
#include <networkit/flow/Dinic.hpp>
namespace NetworKit {
void Dinic::buildResidual() {
    graph->forEdges([&](node u, node v, edgeweight w) {
        residual.addEdge(u, v, w);
        residual.addEdge(v, u, 0.0);
    });
}

bool Dinic::bfs() {
    std::ranges::fill(level, -1);
    std::queue<node> queue;
    level[source] = 0;
    queue.push(source);
    while (!queue.empty()) {
        const node u = queue.front();
        queue.pop();
        for (const auto v : residual.neighborRange(u)) {
            if (residual.weight(u, v) > 0 && level[v] < 0) {
                level[v] = level[u] + 1;
                queue.push(v);
            }
        }
    }
    return level[target] >= 0;
}

edgeweight Dinic::dfs(node u, edgeweight flow) {
    if (u == target || flow == 0)
        return flow;
    const auto neighbor = residual.neighborRange(u).begin();
    while (neighbor != residual.neighborRange(u).end()) {
        const node v = *neighbor;
        const auto remainingCapacity = residual.weight(u, v);
        if (remainingCapacity > 0 && level[v] == level[u] + 1) {
            edgeweight pushed = dfs(v, std::min(flow, remainingCapacity));
            if (pushed > 0.0) {
                residual.setWeight(u, v, remainingCapacity - pushed);
                residual.setWeight(v, u, residual.weight(v, u) + pushed);
                return pushed;
            }
        }
    }
    return 0.0;
}

void Dinic::run() {
    buildResidual();
    edgeweight flow = 0.0;
    level.resize(residual.upperNodeIdBound());
    ptr.resize(residual.upperNodeIdBound());
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
