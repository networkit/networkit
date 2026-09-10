/*
 * TopologicalSortImpl.hpp
 *
 *  Created on: 10.09.2026
 *      Author: NetworKit contributors
 */
#ifndef NETWORKIT_GRAPH_TOPOLOGICAL_SORT_IMPL_HPP_
#define NETWORKIT_GRAPH_TOPOLOGICAL_SORT_IMPL_HPP_

#include <algorithm>
#include <sstream>
#include <stack>
#include <stdexcept>

namespace NetworKit {

template <typename GraphT>
GenericTopologicalSort<GraphT>::GenericTopologicalSort(const GraphT &G)
    : G(G), computedNodeIdMap(computeContinuousNodeIds(G)) {
    checkDirected();
}

template <typename GraphT>
GenericTopologicalSort<GraphT>::GenericTopologicalSort(const GraphT &G,
                                                       const NodeIdMapping &nodeIdMap,
                                                       bool checkMapping)
    : G(G), nodeIdMap(&nodeIdMap) {
    checkDirected();
    if (nodeIdMap.size() != G.numberOfNodes())
        throw std::runtime_error(
            "Node id mapping should contain exactly one entry for every node.");
    else if (checkMapping)
        checkNodeIdMap();
}

template <typename GraphT>
typename GenericTopologicalSort<GraphT>::NodeIdMapping
GenericTopologicalSort<GraphT>::computeContinuousNodeIds(const GraphT &G) {
    NodeIdMapping nodeIdMap;
    nodeIdMap.reserve(G.numberOfNodes());

    index continuousId = 0;
    G.forNodes([&](NodeT u) { nodeIdMap.emplace(u, continuousId++); });

    return nodeIdMap;
}

template <typename GraphT>
void GenericTopologicalSort<GraphT>::checkDirected() {
    if (!G.isDirected())
        throw std::runtime_error("Topological sort is defined for directed graphs only.");
}

template <typename GraphT>
void GenericTopologicalSort<GraphT>::checkNodeIdMap() {
    if (!nodeIdMap)
        return;

    const count numberOfNodes = G.numberOfNodes();
    std::vector<bool> checkTable(numberOfNodes);
    for (const auto &entry : *nodeIdMap) {
        const index mappedNode = entry.second;
        if (mappedNode < numberOfNodes && !checkTable[mappedNode])
            checkTable[mappedNode] = true;
        else
            throw std::runtime_error("Node id mapping is not continuous.");
    }
}

template <typename GraphT>
void GenericTopologicalSort<GraphT>::run() {
    reset();

    std::stack<NodeT> nodeStack;

    G.forNodes([&](NodeT u) {
        const index mappedU = mapNode(u);
        if (topSortMark[mappedU] == NodeMark::PERM)
            return;

        nodeStack.push(u);
        do {
            const NodeT v = nodeStack.top();
            const index mappedV = mapNode(v);

            if (topSortMark[mappedV] != NodeMark::NONE) {
                nodeStack.pop();
                if (topSortMark[mappedV] == NodeMark::TEMP) {
                    topSortMark[mappedV] = NodeMark::PERM;
                    topology[current] = v;
                    current--;
                }
            } else {
                topSortMark[mappedV] = NodeMark::TEMP;
                G.forNeighborsOf(v, [&](NodeT w) {
                    const index mappedW = mapNode(w);

                    if (topSortMark[mappedW] == NodeMark::NONE)
                        nodeStack.push(w);
                    else if (topSortMark[mappedW] == NodeMark::TEMP)
                        throw std::runtime_error("Error: the input graph has cycles.");
                });
            }
        } while (!nodeStack.empty());
    });

    hasRun = true;
}

template <typename GraphT>
index GenericTopologicalSort<GraphT>::mapNode(NodeT u) const {
    if (nodeIdMap) {
        const auto it = nodeIdMap->find(u);
        if (it == nodeIdMap->cend()) {
            std::stringstream errorMsg;
            errorMsg << "Node id mapping does not contain node " << u;
            throw std::runtime_error(errorMsg.str());
        }
        return it->second;
    } else if (computedNodeIdMap.has_value()) {
        return computedNodeIdMap.value().at(u);
    } else {
        return static_cast<index>(u);
    }
}

template <typename GraphT>
void GenericTopologicalSort<GraphT>::reset() {
    const count n = G.numberOfNodes();
    if (n == 0)
        throw std::runtime_error("Graph should contain at least one node.");

    if (n != static_cast<count>(topSortMark.size())) {
        topSortMark.resize(n);
        topology.resize(n);
    }
    std::fill(topSortMark.begin(), topSortMark.end(), NodeMark::NONE);
    std::fill(topology.begin(), topology.end(), NodeT{0});
    current = n - 1;
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_TOPOLOGICAL_SORT_IMPL_HPP_
