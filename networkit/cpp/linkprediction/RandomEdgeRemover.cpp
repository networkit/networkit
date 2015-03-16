/*
 * RandomEdgeRemover.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "RandomEdgeRemover.h"

namespace NetworKit {

RandomEdgeRemover::RandomEdgeRemover(const Graph& G) : G(G) {
}

std::pair<Graph, Graph> RandomEdgeRemover::remove(double percentage) {
  if (percentage < 0 || percentage > 1) {
    throw std::invalid_argument("Given percentage is not in [0, 1].");
  }
  Graph remaining(G);
  Graph removed(G.upperNodeIdBound());
  count numEdgesToRemove = 1.0 * percentage * G.numberOfEdges();
  for (index i = 0; i < numEdgesToRemove; ++i) {
    std::pair<node, node> edgeToRemove = remaining.randomEdge();
    remaining.removeEdge(edgeToRemove.first, edgeToRemove.second);
    removed.addEdge(edgeToRemove.first, edgeToRemove.second);
    INFO("Removed edge (", edgeToRemove.first, ", ", edgeToRemove.second, ").");
  }
  return std::make_pair(remaining, removed);
}

} // namespace NetworKit