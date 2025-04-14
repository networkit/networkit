/*  PrimMST.cpp
 *
 *	Created on: 29.03.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <networkit/graph/PrimMSF.hpp>

#include <queue>
#include <vector>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

void PrimMSF::run() {
    // Unweighted graph
    if (!G->isWeighted()) {
        SpanningForest spanningForest(*G);
        spanningForest.run();
        forest = spanningForest.getForest();
        hasRun = true;
        return;
    }
    // Weighted graph
    forest = GraphTools::copyNodes(*G);
    const index numberOfNodes = G->numberOfNodes();
    std::vector<edgeweight> weights(numberOfNodes, infiniteWeight);
    std::vector<node> parents(numberOfNodes, none);
    std::vector<bool> visited(numberOfNodes, false);
    std::priority_queue<std::pair<edgeweight, node>, std::vector<std::pair<edgeweight, node>>,
                        std::greater<>>
        minHeap;

    G->forNodes([&](node startNode) {
        if (visited[startNode]) {
            return;
        }
        minHeap.emplace(nullWeight, startNode);
        weights[startNode] = nullWeight;
        while (!minHeap.empty()) {
            const auto [currentWeight, currentNode] = minHeap.top();
            minHeap.pop();
            if (visited[currentNode]) {
                continue;
            }
            visited[currentNode] = true;
            G->forNeighborsOf(currentNode, [&](node neighbor, edgeweight neighborWeight) {
                if (!visited[neighbor] && weights[neighbor] > neighborWeight) {
                    weights[neighbor] = neighborWeight;
                    minHeap.emplace(neighborWeight, neighbor);
                    parents[neighbor] = currentNode;
                }
            });
        }
    });
    G->forNodes([&](node node) {
        if (parents[node] != none) {
            forest.addEdge(parents[node], node);
            totalWeight += weights[node];
        }
    });
    hasRun = true;
}

} // namespace NetworKit
