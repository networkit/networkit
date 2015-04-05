/*
 * TotalNeighborsIndex.cpp
 *
 *  Created on: 05.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "TotalNeighborsIndex.h"

namespace NetworKit {

TotalNeighborsIndex::TotalNeighborsIndex() : LinkPredictor() {
}

TotalNeighborsIndex::TotalNeighborsIndex(const Graph& G) : LinkPredictor(G), jaccardIndex(G) {
}

double TotalNeighborsIndex::runImpl(node u, node v) {
  std::vector<node> ns = jaccardIndex.getNeighborsUnion(u, v);
  return jaccardIndex.getNeighborsUnion(u, v).size();
}

void TotalNeighborsIndex::setGraph(const Graph& newGraph) {
  jaccardIndex.setGraph(newGraph);
  G = &newGraph;
  validCache = false;
}

} // namespace NetworKit
