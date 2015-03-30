/*
 * RandomEdgePartitioner.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "RandomEdgePartitioner.h"

namespace NetworKit {

RandomEdgePartitioner::RandomEdgePartitioner(const Graph& G) : G(G) {
}

std::pair<Graph, Graph> RandomEdgePartitioner::partitionByPercentage(double percentage) const {
  if (percentage < 0 || percentage > 1) {
    throw std::invalid_argument("Given percentage is not in [0, 1].");
  }
  return partitionByCount(percentage * G.numberOfEdges());
}

std::pair<Graph, Graph> RandomEdgePartitioner::partitionByCount(count numEdges) const {
  if (numEdges > G.numberOfEdges()) {
    throw std::invalid_argument("numEdges > G.numberOfEdges().");
  }
  Graph remaining(G);
  Graph removed(G.upperNodeIdBound());
  for (index i = 0; i < numEdges; ++i) {
    std::pair<node, node> edgeToRemove = remaining.randomEdge();
    remaining.removeEdge(edgeToRemove.first, edgeToRemove.second);
    removed.addEdge(edgeToRemove.first, edgeToRemove.second);
  }
  return std::make_pair(remaining, removed);
}

} // namespace NetworKit