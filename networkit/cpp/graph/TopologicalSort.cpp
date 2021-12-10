/*
 * TopologicalSort.cpp
 *
 *  Created on: 22.11.2021
 *      Author: Fabian Brandt-Tumescheit
 */
#include <networkit/graph/TopologicalSort.hpp>

namespace NetworKit {

TopologicalSort::TopologicalSort(const Graph &G) : G(&G) {
    if (!G.isDirected())
        throw std::runtime_error("Topological sort is defined for directed graphs only.");
}

void TopologicalSort::run() {
    reset();

    std::stack<node> nodeStack;

    G->forNodes([&](node u) {
        if (topSortMark[u] == NodeMark::PERM)
            return;

        nodeStack.push(u);
        do {
            node v = nodeStack.top();
            if (topSortMark[v] != NodeMark::NONE) {
                nodeStack.pop();
                if (topSortMark[v] == NodeMark::TEMP) {
                    topSortMark[v] = NodeMark::PERM;
                    topology[current] = v;
                    current--;
                }
            } else {
                topSortMark[v] = NodeMark::TEMP;
                G->forNeighborsOf(v, [&](node w) {
                    if (topSortMark[w] == NodeMark::NONE)
                        nodeStack.push(w);
                    else if (topSortMark[w] == NodeMark::TEMP)
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