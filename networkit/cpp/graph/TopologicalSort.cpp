/*
 * TopologicalSort.cpp
 *
 *  Created on: 22.11.2021
 *      Author: Fabian Brandt-Tumescheit
 */
#include <stdexcept>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/TopologicalSort.hpp>

namespace NetworKit {

TopologicalSort::TopologicalSort(const Graph &G)
    : G(G), computedNodeIdMap(GraphTools::getContinuousNodeIds(G)) {
    checkDirected();
}

TopologicalSort::TopologicalSort(const Graph &G, const std::unordered_map<node, node> &nodeIdMap,
                                 bool checkMapping)
    : G(G), nodeIdMap(&nodeIdMap) {
    checkDirected();
    if (nodeIdMap.size() != G.numberOfNodes())
        throw std::runtime_error(
            "Node id mapping should contain exactly one entry for every node.");
    else if (checkMapping)
        checkNodeIdMap();
}

void TopologicalSort::checkDirected() {
    if (!G.isDirected())
        throw std::runtime_error("Topological sort is defined for directed graphs only.");
}

void TopologicalSort::checkNodeIdMap() {
    if (!nodeIdMap)
        return;

    size_t numberOfNodes = G.numberOfNodes();
    std::vector<bool> checkTable(numberOfNodes);
    for (auto &[origNode, mappedNode] : *nodeIdMap) {
        if (mappedNode < numberOfNodes && !checkTable[mappedNode])
            checkTable[mappedNode] = true;
        else
            throw std::runtime_error("Node id mapping is not continuous.");
    }
}

void TopologicalSort::run() {
    reset();

    std::stack<node> nodeStack;

    G.forNodes([&](node u) {
        node mappedU = mapNode(u);
        if (topSortMark[mappedU] == NodeMark::PERM)
            return;

        nodeStack.push(u);
        do {
            node v = nodeStack.top();
            node mappedV = mapNode(v);

            if (topSortMark[mappedV] != NodeMark::NONE) {
                nodeStack.pop();
                if (topSortMark[mappedV] == NodeMark::TEMP) {
                    topSortMark[mappedV] = NodeMark::PERM;
                    topology[current] = v;
                    current--;
                }
            } else {
                topSortMark[mappedV] = NodeMark::TEMP;
                G.forNeighborsOf(v, [&](node w) {
                    node mappedW = mapNode(w);

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

node TopologicalSort::mapNode(node u) {
    if (nodeIdMap) {
        auto it = nodeIdMap->find(u);
        if (it == nodeIdMap->cend()) {
            std::stringstream errorMsg;
            errorMsg << "Node id mapping does not contain node " << u;
            throw std::runtime_error(errorMsg.str());
        }
        return it->second;
    } else if (computedNodeIdMap.has_value())
        return computedNodeIdMap.value().at(u);
    else
        return u;
}

void TopologicalSort::reset() {
    // Reset node marks
    count n = G.numberOfNodes();
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
