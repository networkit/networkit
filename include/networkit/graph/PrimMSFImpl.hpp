/*  PrimMSFImpl.hpp
 *
 *	Created on: 11.09.2026
 *  Authors: NetworKit contributors
 *
 */
#ifndef NETWORKIT_GRAPH_PRIM_MSF_IMPL_HPP_
#define NETWORKIT_GRAPH_PRIM_MSF_IMPL_HPP_

#include <cassert>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include <tlx/container/d_ary_heap.hpp>

namespace NetworKit {

template <typename GraphT>
void GenericPrimMSF<GraphT>::run() {
    using HeapElement = std::pair<EdgeWeightT, NodeT>;

    // Unweighted graph
    if (!this->G->isWeighted()) {
        GenericSpanningForest<GraphT> spanningForest(*this->G);
        spanningForest.run();
        this->forest = spanningForest.getForest();
        this->hasRun = true;
        return;
    }
    // Weighted graph
    this->forest = GenericSpanningForest<GraphT>::copyNodes(*this->G);
    const index nodeIdBound = static_cast<index>(this->G->upperNodeIdBound());
    std::vector<EdgeWeightT> weights(nodeIdBound, infiniteWeight);
    std::vector<NodeT> parents(nodeIdBound, NullNodeId<NodeT>);
    std::vector<bool> visited(nodeIdBound, false);
    tlx::d_ary_heap<HeapElement, 2, std::less<HeapElement>> minHeap;

    const auto nodeIndex = [](NodeT u) {
        if constexpr (std::is_signed_v<NodeT>) {
            assert(u >= 0);
        }
        return static_cast<index>(u);
    };

    this->G->forNodes([&](NodeT startNode) {
        const index startNodeIndex = nodeIndex(startNode);
        if (visited[startNodeIndex]) {
            return;
        }
        minHeap.push({EdgeWeightT{0}, startNode});
        weights[startNodeIndex] = EdgeWeightT{0};
        parents[startNodeIndex] = startNode;
        while (!minHeap.empty()) {
            const auto pair = minHeap.top();
            const EdgeWeightT currentWeight = pair.first;
            const NodeT currentNode = pair.second;
            const index currentNodeIndex = nodeIndex(currentNode);
            minHeap.pop();
            if (visited[currentNodeIndex]) {
                continue;
            }
            if (const NodeT parentNode = parents[currentNodeIndex]; currentNode != parentNode) {
                this->forest.addEdge(parentNode, currentNode);
                totalWeight += static_cast<edgeweight>(currentWeight);
            }
            visited[currentNodeIndex] = true;
            for (const auto &[neighbor, neighborWeight] :
                 this->G->weightNeighborRange(currentNode)) {
                const index neighborIndex = nodeIndex(neighbor);
                if (!visited[neighborIndex] && weights[neighborIndex] > neighborWeight) {
                    weights[neighborIndex] = neighborWeight;
                    minHeap.push({neighborWeight, neighbor});
                    parents[neighborIndex] = currentNode;
                }
            }
        }
    });
    this->hasRun = true;
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_PRIM_MSF_IMPL_HPP_
