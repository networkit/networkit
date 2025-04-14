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
    std::vector<bool> visited(numberOfNodes, false);
    std::priority_queue<std::pair<edgeweight, node>, std::vector<std::pair<edgeweight, node>>, std::greater<>> minHeap;

    G->forNodes([&](node startVertex) {
        if (visited[startVertex]) {
            return;
        }
        minHeap.emplace(nullWeight, startVertex);
        weights[startVertex] = nullWeight;
        while(!minHeap.empty()) {
            auto [currentWeight, currentNode] = minHeap.top();
            minHeap.pop();
            if(visited[currentNode]) {
                continue;
            }
            visited[currentNode] = true;
            G->forNeighborsOf( currentNode, [&](node neighbor, edgeweight neighborWeight) {
                if(!visited[neighbor] && weights[neighbor] > neighborWeight) {
                    weights[neighbor] = neighborWeight;
                    minHeap.emplace(neighborWeight, neighbor);
                    forest.addEdge(currentNode, neighbor);
                    totalWeight += neighborWeight;
                }
            });
        }
    });
    hasRun = true;
}


}
