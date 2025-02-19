/*  FloydWarshall.cpp
 *
 *	Created on: 15.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <networkit/distance/FloydWarshall.hpp>

namespace NetworKit {

FloydWarshall::FloydWarshall(const Graph &G, double densityThreshold, node maximumNumberOfNodes)
    : graph(&G), numberOfNodes(G.numberOfNodes()) {

    if (!G.isWeighted()) {
        throw std::runtime_error("The input graph is unweighted!");
    }
    if (densityThreshold <= 0.0 || densityThreshold > 1.0) {
        throw std::invalid_argument("Invalid density threshold. Must be in range (0, 1].");
    }
    if (maximumNumberOfNodes < 1) {
        throw std::invalid_argument("Invalid maximum node count. Must be at least 1.");
    }

    const index maximumNumberOfEdges = G.isDirected() ? maximumNumberOfNodes*(maximumNumberOfNodes-1) : maximumNumberOfNodes*(maximumNumberOfNodes-1)/2;
    const double density = static_cast<double>(G.numberOfEdges()) / maximumNumberOfEdges;

    if (density < densityThreshold) {
        throw std::domain_error("Graph density is below user-defined density-threshold of: " + std::to_string(densityThreshold));
    }
    if (numberOfNodes > maximumNumberOfNodes) {
        throw std::domain_error("Graph size exceeds user-defined max-node-limit of: " + std::to_string(maximumNumberOfNodes));
    }
}

void FloydWarshall::tagNegativeCycles() {
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

std::vector<std::vector<edgeweight>> FloydWarshall::getAllDistances() const {
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
