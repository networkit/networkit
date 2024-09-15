/*
 * TopologicalSort.cpp
 *
 *  Created on: 22.11.2021
 *      Author: Fabian Brandt-Tumescheit
 */
#include <stdexcept>
#include "networkit/graph/GraphTools.hpp"
#include <networkit/graph/TopologicalSort.hpp>

namespace NetworKit {

TopologicalSort::TopologicalSort(const Graph &G)
    : G(&G), computedNodeIdMap(GraphTools::getContinuousNodeIds(G)) {
    checkDirected();
    if (G.upperNodeIdBound() != G.numberOfNodes() - 1) {
        computedNodeIdMap = GraphTools::getContinuousNodeIds(G);
        nodeIdMap = &computedNodeIdMap.value();
    }
}

TopologicalSort::TopologicalSort(const Graph &G, std::unordered_map<node, node> &nodeIdMap,
                                 bool checkMapping)
    : G(&G), computedNodeIdMap({}), nodeIdMap(&nodeIdMap) {
    checkDirected();
    size_t numberOfNodes = G.numberOfNodes();
    if (nodeIdMap.size() != numberOfNodes)
        throw std::runtime_error(
            "Node id mapping should contain exactly one entry for every node.");

    if (checkMapping)
        checkNodeIdMap();
}

void TopologicalSort::checkDirected() {
    if (!G->isDirected())
        throw std::runtime_error("Topological sort is defined for directed graphs only.");
}

void TopologicalSort::checkNodeIdMap() {
    size_t numberOfNodes = G->numberOfNodes();
    std::vector<bool> checkTable(numberOfNodes);
    for (auto it = nodeIdMap->begin(); it != nodeIdMap->end(); it++) {
        node mappedNode = it->second;
        if (mappedNode < numberOfNodes && !checkTable[mappedNode])
            checkTable[mappedNode] = true;
        else
            throw std::runtime_error("Node id mapping is not continuous.");
    }
}

void TopologicalSort::run() {
    reset();

    std::stack<node> nodeStack;

    G->forNodes([&](node u) {
        node mappedU;
        if (nodeIdMap != nullptr)
            mappedU = (*nodeIdMap)[u];
        else
            mappedU = u;

        if (topSortMark[mappedU] == NodeMark::PERM)
            return;

        nodeStack.push(u);
        do {
            node v = nodeStack.top();
            node mappedV;
            if (nodeIdMap != nullptr)
                mappedV = (*nodeIdMap)[v];
            else
                mappedV = v;

            if (topSortMark[mappedV] != NodeMark::NONE) {
                nodeStack.pop();
                if (topSortMark[mappedV] == NodeMark::TEMP) {
                    topSortMark[mappedV] = NodeMark::PERM;
                    topology[current] = v;
                    current--;
                }
            } else {
                topSortMark[mappedV] = NodeMark::TEMP;
                G->forNeighborsOf(v, [&](node w) {
                    node mappedW;
                    if (nodeIdMap != nullptr)
                        mappedW = (*nodeIdMap)[w];
                    else
                        mappedW = w;

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

void TopologicalSort::reset() {
    // Reset node marks
    count n = G->numberOfNodes();
    if (n == 0)
        throw std::runtime_error("Graph should contain at least one node.");

    if (n != static_cast<count>(topSortMark.size())) {
        topSortMark.resize(n);
        topology.resize(n);
    }
    std::fill(topSortMark.begin(), topSortMark.end(), NodeMark::NONE);
    std::fill(topology.begin(), topology.end(), 0);
    // Reset current
    current = n - 1;
}

} // namespace NetworKit
