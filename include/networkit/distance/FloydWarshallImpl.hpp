/*  FloydWarshallImpl.hpp
 *
 *	Created on: 15.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_DISTANCE_FLOYD_WARSHALL_IMPL_HPP_
#define NETWORKIT_DISTANCE_FLOYD_WARSHALL_IMPL_HPP_

#include <limits>
#include <stdexcept>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/NumericTools.hpp>

namespace NetworKit {

template <typename GraphT>
GenericFloydWarshall<GraphT>::GenericFloydWarshall(const GraphT &G) : graph(&G) {
    if (!G.isWeighted()) {
        throw std::runtime_error("The input graph is unweighted!");
    }
}

template <typename GraphT>
void GenericFloydWarshall<GraphT>::tagNegativeCycles() {
    graph->forNodes([&](NodeT w) {
        if (distances[w][w] >= 0.0)
            return;
        nodesInNegativeCycle[w] = 1;
        graph->forNodes([&](NodeT u) {
            if (distances[u][w] == infiniteDistance)
                return;
            graph->forNodes([&](NodeT v) {
                if (distances[w][v] != infiniteDistance) {
                    nodesInNegativeCycle[u] = 1;
                    nodesInNegativeCycle[v] = 1;
                    distances[u][v] = -std::numeric_limits<DistanceT>::infinity();
                    pathMatrix[u][v] = NullNodeId<NodeT>;
                }
            });
        });
    });
}

template <typename GraphT>
void GenericFloydWarshall<GraphT>::run() {
    const index numberOfNodes = graph->numberOfNodes();
    distances.resize(numberOfNodes, std::vector<DistanceT>(numberOfNodes, infiniteDistance));
    nodesInNegativeCycle.resize(numberOfNodes);
    pathMatrix.resize(numberOfNodes, std::vector(numberOfNodes, NullNodeId<NodeT>));
    hops.resize(numberOfNodes, std::vector(numberOfNodes, none));

    graph->forNodes([&](NodeT u) {
        distances[u][u] = 0.0;
        pathMatrix[u][u] = u;
        hops[u][u] = 0;
    });

    graph->forNodes([&](NodeT u) {
        for (const NodeT v : graph->neighborRange(u)) {
            distances[u][v] = static_cast<DistanceT>(graph->weight(u, v));
            pathMatrix[u][v] = v;
            hops[u][v] = 1;
        }
    });

    graph->forNodes([&](NodeT intermediate) {
        graph->parallelForNodes([&](NodeT source) {
            if (distances[source][intermediate] == infiniteDistance)
                return;
            graph->forNodes([&](NodeT target) {
                if (distances[intermediate][target] == infiniteDistance) {
                    return;
                }
                const DistanceT candidateDistance =
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

template <typename GraphT>
typename GenericFloydWarshall<GraphT>::DistanceT
GenericFloydWarshall<GraphT>::getDistance(NodeT source, NodeT target) const {
    assureFinished();
    return distances[source][target];
}

template <typename GraphT>
bool GenericFloydWarshall<GraphT>::isNodeInNegativeCycle(NodeT u) const {
    assureFinished();
    return nodesInNegativeCycle[u] == 1;
}

template <typename GraphT>
std::vector<typename GenericFloydWarshall<GraphT>::NodeT>
GenericFloydWarshall<GraphT>::getNodesOnShortestPath(NodeT source, NodeT target) const {
    assureFinished();
    if (pathMatrix[source][target] == NullNodeId<NodeT>) {
        return {};
    }
    std::vector<NodeT> path;
    NodeT currentNode = source;

    while (currentNode != target) {
        if (currentNode == NullNodeId<NodeT>)
            return {};
        path.push_back(currentNode);
        currentNode = pathMatrix[currentNode][target];
    }
    path.push_back(target);
    return path;
}

template <typename GraphT>
const std::vector<std::vector<typename GenericFloydWarshall<GraphT>::DistanceT>> &
GenericFloydWarshall<GraphT>::getDistances() const {
    assureFinished();
    return distances;
}

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_FLOYD_WARSHALL_IMPL_HPP_
