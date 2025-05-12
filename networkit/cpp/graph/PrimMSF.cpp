/*  PrimMST.cpp
 *
 *	Created on: 29.03.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <networkit/graph/PrimMSF.hpp>

#include <vector>
#include <tlx/container/d_ary_heap.hpp>
#include <networkit/auxiliary/VectorComparator.hpp>
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
    tlx::d_ary_heap<std::pair<edgeweight, node>, 2, std::less<std::pair<edgeweight, node>>> minHeap;

    G->forNodes([&](node startNode) {
        if (visited[startNode]) {
            return;
        }
        minHeap.push({nullWeight, startNode});
        weights[startNode] = nullWeight;
        parents[startNode] = startNode;
        while (!minHeap.empty()) {
            const auto pair = minHeap.top();
            const edgeweight currentWeight = pair.first;
            const node currentNode = pair.second;
            minHeap.pop();
            if (visited[currentNode]) {
                continue;
            }
            if (const node parentNode = parents[currentNode]; currentNode != parentNode) {
                forest.addEdge(parentNode, currentNode);
                totalWeight += currentWeight;
            }
            visited[currentNode] = true;
            G->forNeighborsOf(currentNode, [&](node neighbor, edgeweight neighborWeight) {
                if (!visited[neighbor] && weights[neighbor] > neighborWeight) {
                    weights[neighbor] = neighborWeight;
                    minHeap.push({neighborWeight, neighbor});
                    parents[neighbor] = currentNode;
                }
            });
        }
    });
    hasRun = true;
}

} // namespace NetworKit
