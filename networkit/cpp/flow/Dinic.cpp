/*  Dinic.cpp
 *
 *	Created on: 20.06.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <networkit/flow/Dinic.hpp>

namespace NetworKit {
void Dinic::buildResidual() {
    graph->forEdges([&](node u, node v, edgeweight w) {
        residual.addEdge(u, v, w);
        residual.addEdge(v, u, 0.0);
    });
}

bool Dinic::bfs() {
    std::ranges::fill(level.begin(), level.end(), -1);
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

void Dinic::run() {
    buildResidual();
    edgeweight flow = 0.0;
    const edgeweight INF = std::numeric_limits<edgeweight>::max();
    level.resize(residual.upperNodeIdBound());
    ptr.resize(residual.upperNodeIdBound());
    // Main loop: while there is an augmenting path
    while (bfs()) {
        std::ranges::fill(ptr.begin(), ptr.end(), 0);
        // Send blocking flow
        while (auto pushed = dfs(source, std::numeric_limits<edgeweight>::max())) {
            flow += pushed;
        }
    }
    maxFlow = flow;
    hasRun = true;
}
} // namespace NetworKit
