/*
 * CommonNeighborsIndex.cpp
 *
 *  Created on: 06.12.2014
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "CommonNeighborsIndex.h"

namespace NetworKit {

CommonNeighborsIndex::CommonNeighborsIndex() {
}

CommonNeighborsIndex::CommonNeighborsIndex(const Graph& G) : LinkPredictor(G) {
}

double CommonNeighborsIndex::runImpl(node u, node v) {
  return getCommonNeighbors(u, v).size();
}

std::vector<node> CommonNeighborsIndex::getCommonNeighbors(node u, node v) const {
  std::vector<node> uNeighbors = G->neighbors(u);
  std::vector<node> vNeighbors = G->neighbors(v);
  std::vector<node> commonNeighbors;
  // We have no guarantee that the neighbor-vectors are sorted so we have to
  // sort them in order for set_intersection to work properly.
  std::sort(uNeighbors.begin(), uNeighbors.end());
  std::sort(vNeighbors.begin(), vNeighbors.end());
  std::set_intersection(uNeighbors.begin(), uNeighbors.end(), vNeighbors.begin(),
    vNeighbors.end(), std::back_inserter(commonNeighbors));
  return commonNeighbors;
}

} // namespace NetworKit