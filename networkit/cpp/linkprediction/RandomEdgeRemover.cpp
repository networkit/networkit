/*
 * RandomEdgeRemover.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "../auxiliary/Log.h"

#include "RandomEdgeRemover.h"

namespace NetworKit {

RandomEdgeRemover::RandomEdgeRemover(const Graph& G) : G(G) {
}

std::pair<Graph, Graph> RandomEdgeRemover::remove(double percentage) {
  if (percentage < 0 || percentage > 1) {
    throw std::invalid_argument("Given percentage is not in [0, 1].");
  }
  // TODO: Check if this is a deep copy
  Graph remaining(G);
  INFO("[RandomEdgeRemover] Given Graph has ", remaining.numberOfEdges(), " edges.\n");
  Graph removed(G.upperNodeIdBound());
  count numEdgesToRemove = 1.0 * percentage * G.numberOfEdges();
  //std::vector<std::pair<node, node>> edgesToRemove = remaining.randomEdges(numEdgesToRemove);
  INFO("[RandomEdgeRemover] There are ", numEdgesToRemove, " edges to remove.\n");
  for (index i = 0; i < numEdgesToRemove; ++i) {
    std::pair<node, node> edgeToRemove = remaining.randomEdge();
    INFO("[RandomEdgeRemover] Removing edge (", edgeToRemove.first, ", ", edgeToRemove.second, ").\n");
    remaining.removeEdge(edgeToRemove.first, edgeToRemove.second);
    removed.addEdge(edgeToRemove.first, edgeToRemove.second);
  }
  INFO("G: ",)
  INFO("")
  return std::make_pair(remaining, removed);
}

} // namespace NetworKit