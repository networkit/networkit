/*  FloydWarshall.cpp
 *
 *	Created on: 15.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <networkit/distance/FloydWarshall.hpp>

namespace NetworKit {

FloydWarshall::FloydWarshall(const Graph &G)
    : graph(&G) {

    if (!G.isWeighted()) {
        throw std::runtime_error("The input graph is unweighted!");
    }
}

void FloydWarshall::tagNegativeCycles() {
    const index numberOfNodes = graph->numberOfNodes();
    for (node w = 0; w < numberOfNodes; ++w) {
        if (!(distances[w][w] < 0.0))
            continue;
        nodesInNegativeCycle[w] = 1;
        for (node u = 0; u < numberOfNodes; ++u) {
            if (distances[u][w] == std::numeric_limits<edgeweight>::max())
                continue;
            for (node v = 0; v < numberOfNodes; ++v) {
                if (distances[w][v] != std::numeric_limits<edgeweight>::max()) {
                    nodesInNegativeCycle[u] = 1;
                    nodesInNegativeCycle[v] = 1;
                    distances[u][v] = -std::numeric_limits<edgeweight>::infinity();
                    pathMatrix[u][v] = none;
                }
            }
        }
    }
}

void FloydWarshall::run() {
    const index numberOfNodes = graph->numberOfNodes();
    distances = std::vector(numberOfNodes,
                            std::vector(numberOfNodes, std::numeric_limits<edgeweight>::max()));
    nodesInNegativeCycle = std::vector<uint8_t>(numberOfNodes);
    pathMatrix = std::vector(numberOfNodes, std::vector(numberOfNodes, none));

    for (node u = 0; u < numberOfNodes; ++u) {
        distances[u][u] = 0.0;
        pathMatrix[u][u] = u;
    }

    for (node u = 0; u < numberOfNodes; ++u) {
        for (const node v : graph->neighborRange(u)) {
            distances[u][v] = graph->weight(u, v);
            pathMatrix[u][v] = v;
        }
    }

    for (node w = 0; w < numberOfNodes; ++w) {
        for (node u = 0; u < numberOfNodes; ++u) {
            if (distances[u][w] == std::numeric_limits<edgeweight>::max())
                continue;
            for (node v = 0; v < numberOfNodes; ++v) {
                if (distances[w][v] != std::numeric_limits<edgeweight>::max()
                    && distances[u][w] + distances[w][v] < distances[u][v]) {
                    distances[u][v] = distances[u][w] + distances[w][v];
                    pathMatrix[u][v] = pathMatrix[u][w];
                }
            }
        }
    }

    tagNegativeCycles();
    hasRun = true;
}

edgeweight FloydWarshall::getDistance(node source, node target) const {
    assureFinished();
    return distances[source][target];
}

const std::vector<std::vector<edgeweight>> & FloydWarshall::getAllDistances() const & {
    assureFinished();
    return distances;
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
