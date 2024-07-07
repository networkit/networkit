/*
 * TopologicalSort.cpp
 *
 *  Created on: 22.11.2021
 *      Author: Fabian Brandt-Tumescheit
 */
#include <stdexcept>
#include <networkit/graph/TopologicalSort.hpp>

namespace NetworKit {

TopologicalSort::TopologicalSort(const Graph &G) : G(&G), nodeIdMap(nullptr) {
    checkDirected();
}

TopologicalSort::TopologicalSort(const Graph &G,
                                 const std::unordered_map<node, node> &nodeIdMap)
    : G(&G), nodeIdMap(&nodeIdMap) {
    checkDirected();
}

void TopologicalSort::checkDirected() {
    if (!G->isDirected())
        throw std::runtime_error("Topological sort is defined for directed graphs only.");
}

void TopologicalSort::run() {
    reset();

    std::stack<node> nodeStack;

    try {
        G->forNodes([&](node u) {
            node mappedU;
            if (nodeIdMap != nullptr)
                mappedU = nodeIdMap->at(u);
            else
                mappedU = u;

            if (topSortMark.at(mappedU) == NodeMark::PERM)
                return;

            nodeStack.push(u);
            do {
                node v = nodeStack.top();
                node mappedV;
                if (nodeIdMap != nullptr)
                    mappedV = nodeIdMap->at(v);
                else
                    mappedV = v;

                if (topSortMark.at(mappedV) != NodeMark::NONE) {
                    nodeStack.pop();
                    if (topSortMark.at(mappedV) == NodeMark::TEMP) {
                        topSortMark.at(mappedV) = NodeMark::PERM;
                        topology[current] = v;
                        current--;
                    }
                } else {
                    topSortMark.at(mappedV) = NodeMark::TEMP;
                    G->forNeighborsOf(v, [&](node w) {
                        node mappedW;
                        if (nodeIdMap != nullptr)
                            mappedW = nodeIdMap->at(w);
                        else
                            mappedW = w;

                        if (topSortMark.at(mappedW) == NodeMark::NONE)
                            nodeStack.push(w);
                        else if (topSortMark.at(mappedW) == NodeMark::TEMP)
                            throw std::runtime_error("Error: the input graph has cycles.");
                    });
                }
            } while (!nodeStack.empty());
        });

        hasRun = true;
    } catch (const std::out_of_range &) {
        if (nodeIdMap != nullptr) {
            throw std::runtime_error("Error: node id mapping does not contain all nodes");
        }

        throw std::runtime_error("Error: node ids are not continuous");
    }
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
