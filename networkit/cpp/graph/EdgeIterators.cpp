/*
 * Graph.cpp
 *
 *  Created on: 16.09.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#include <networkit/graph/EdgeIterators.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

template <>
bool EdgeIteratorBase<Graph>::validEdge() const noexcept {
    return G->isDirected() || (*nodeIter <= G->getIthNeighbor(Unsafe{}, *nodeIter, i));
}

template <>
Edge EdgeTypeIterator<Graph, Edge>::operator*() const noexcept {
    assert(nodeIter != G->nodeRange().end());
    return Edge(*nodeIter, G->getIthNeighbor(Unsafe{}, *nodeIter, i));
}

template <>
WeightedEdge EdgeTypeIterator<Graph, WeightedEdge>::operator*() const noexcept {
    assert(nodeIter != G->nodeRange().end());
    return WeightedEdge(*nodeIter, G->getIthNeighbor(Unsafe{}, *nodeIter, i),
                        G->getIthNeighborWeight(Unsafe{}, *nodeIter, i));
}

} // namespace NetworKit
