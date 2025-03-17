/*  FloydWarshall.cpp
 *
 *	Created on: 15.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/distance/FloydWarshall.hpp>

namespace NetworKit {

FloydWarshall::FloydWarshall(const Graph &G) : graph(&G) {
    if (!G.isWeighted()) {
        throw std::runtime_error("The input graph is unweighted!");
    }
}

void FloydWarshall::tagNegativeCycles() {
    graph->forNodes([&](node w) {
        if (distances[w][w] >= 0.0)
            return;
        nodesInNegativeCycle[w] = 1;
        graph->forNodes([&](node u) {
            if (distances[u][w] == infiniteDistance)
                return;
            graph->forNodes([&](node v) {
                if (distances[w][v] != infiniteDistance) {
                    nodesInNegativeCycle[u] = 1;
                    nodesInNegativeCycle[v] = 1;
                    distances[u][v] = -std::numeric_limits<edgeweight>::infinity();
                    pathMatrix[u][v] = none;
                }
            });
        });
    });
}

void FloydWarshall::run() {
    const index numberOfNodes = graph->numberOfNodes();
    distances.resize(numberOfNodes, std::vector<edgeweight>(numberOfNodes, infiniteDistance));
    nodesInNegativeCycle.resize(numberOfNodes);
    pathMatrix.resize(numberOfNodes, std::vector(numberOfNodes, none));
    hops.resize(numberOfNodes, std::vector(numberOfNodes, none));

    graph->forNodes([&](node u) {
        distances[u][u] = 0.0;
        pathMatrix[u][u] = u;
        hops[u][u] = 0;
    });

    graph->forNodes([&](node u) {
        for (const node v : graph->neighborRange(u)) {
            distances[u][v] = graph->weight(u, v);
            pathMatrix[u][v] = v;
            hops[u][v] = 1;
        }
    });

    graph->forNodes([&](node intermediate) {
        graph->parallelForNodes([&](node source) {
            if (distances[source][intermediate] == infiniteDistance)
                return;
            graph->forNodes([&](node target) {
                if (distances[intermediate][target] == infiniteDistance) {
                    return;
                }
                const edgeweight candidateDistance =
                    distances[source][intermediate] + distances[intermediate][target];
                const count candidateHops = hops[source][intermediate] + hops[intermediate][target];
                if (candidateDistance < distances[source][target]) {
                    distances[source][target] = candidateDistance;
                    hops[source][target] = candidateHops;
                    pathMatrix[source][target] = pathMatrix[source][intermediate];
                }
                if (Aux::NumericTools::equal(candidateDistance, distances[source][target])
                    && candidateHops < hops[source][target]) {
                    hops[source][target] = candidateHops;
                    pathMatrix[source][target] = pathMatrix[source][intermediate];
                }
            });
        });
    });

    tagNegativeCycles();
    hasRun = true;
}

edgeweight FloydWarshall::getDistance(node source, node target) const {
    assureFinished();
    return distances[source][target];
}

bool FloydWarshall::isNodeInNegativeCycle(node u) const {
    assureFinished();
    return nodesInNegativeCycle[u] == 1;
}

std::vector<node> FloydWarshall::getNodesOnShortestPath(node source, node target) const {
    assureFinished();
    if (pathMatrix[source][target] == none) {
        return {};
    }
    std::vector<node> path;
    node currentNode = source;

    while (currentNode != target) {
        if (currentNode == none)
            return {};
        path.push_back(currentNode);
        currentNode = pathMatrix[currentNode][target];
    }
    path.push_back(target);
    return path;
}

} // namespace NetworKit
