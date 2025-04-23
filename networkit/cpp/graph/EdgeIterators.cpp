/*
 * Graph.cpp
 *
 *  Created on: 16.09.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#include <networkit/graph/DHBGraph.hpp>
#include <networkit/graph/EdgeIterators.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

template <>
bool EdgeIteratorBase<Graph>::validEdge() const noexcept {
    return G->isDirected() || (*nodeIter <= G->getIthNeighbor(Unsafe{}, *nodeIter, i));
}

template <>
bool EdgeIteratorBase<DHBGraph>::validEdge() const noexcept {
    return G->isDirected() || (*nodeIter <= G->getIthNeighbor(*nodeIter, i));
}

template <>
Edge EdgeTypeIterator<Graph, Edge>::operator*() const noexcept {
    assert(nodeIter != G->nodeRange().end());
    return Edge(*nodeIter, G->getIthNeighbor(Unsafe{}, *nodeIter, i));
}

template <>
Edge EdgeTypeIterator<DHBGraph, Edge>::operator*() const noexcept {
    assert(nodeIter != G->nodeRange().end());
    return Edge(*nodeIter, G->getIthNeighbor(*nodeIter, i));
}

template <>
WeightedEdge EdgeTypeIterator<Graph, WeightedEdge>::operator*() const noexcept {
    assert(nodeIter != G->nodeRange().end());
    return WeightedEdge(*nodeIter, G->getIthNeighbor(Unsafe{}, *nodeIter, i),
                        G->getIthNeighborWeight(Unsafe{}, *nodeIter, i));
}

template <>
WeightedEdge EdgeTypeIterator<DHBGraph, WeightedEdge>::operator*() const noexcept {
    assert(nodeIter != G->nodeRange().end());
    return WeightedEdge(*nodeIter, G->getIthNeighbor(*nodeIter, i),
                        G->isWeighted() ? G->getIthNeighborWeight(*nodeIter, i) : 1);
}

} // namespace NetworKit
