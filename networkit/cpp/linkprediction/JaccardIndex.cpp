/*
 * JaccardIndex.cpp
 *
 *  Created on: 23.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "JaccardIndex.h"

namespace NetworKit {

JaccardIndex::JaccardIndex() : commonNeighborsIndex() {
}

JaccardIndex::JaccardIndex(const Graph& G) : LinkPredictor(G), commonNeighborsIndex() {
}

double JaccardIndex::runImpl(node u, node v) {
  commonNeighborsIndex.setGraph(*G);
  return commonNeighborsIndex.run(u, v) / getNeighborsUnion(u, v).size();
}

std::vector<node> JaccardIndex::getNeighborsUnion(node u, node v) const {
  if (G == nullptr) {
    throw std::logic_error("Set a graph first.");
  } else if (!G->hasNode(u) || !G->hasNode(v)) {
    throw std::invalid_argument("Invalid node provided.");
  }
  std::vector<node> uNeighbors = G->neighbors(u);
  std::vector<node> vNeighbors = G->neighbors(v);
  std::vector<node> neighborsUnion;
  // We have no guarantee that the neighbor-vectors are sorted so we have to
  // sort them in order for set_intersection to work properly.
  std::sort(uNeighbors.begin(), uNeighbors.end());
  std::sort(vNeighbors.begin(), vNeighbors.end());
  std::set_union(uNeighbors.begin(), uNeighbors.end(), vNeighbors.begin(),
    vNeighbors.end(), std::back_inserter(neighborsUnion));
  return neighborsUnion;
}


} // namespace NetworKit
