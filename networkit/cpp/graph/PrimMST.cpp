/*  PrimMST.cpp
 *
 *	Created on: 29.03.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <networkit/graph/PrimMST.hpp>

#include <vector>
#include <queue>

namespace NetworKit {
void PrimMST::run() {
    const index numberOfNodes = G->numberOfNodes();
    std::vector<edgeweight> weights(numberOfNodes, infiniteWeight);
    std::vector<index> parents(numberOfNodes, none);
    std::vector<bool> visited(numberOfNodes, false);
    std::priority_queue<std::pair<edgeweight, node>, std::vector<std::pair<edgeweight, node>, std::greater<>>> min_heap;
    node start_vertex =0;
    while(!G->hasNode(start_vertex)) {
        ++start_vertex;
    }
    min_heap.push({nullWeight, start_vertex});
    weights[start_vertex] = nullWeight;


}


}
