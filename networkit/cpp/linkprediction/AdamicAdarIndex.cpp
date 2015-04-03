/*
 * AdamicAdarIndex.cpp
 *
 *  Created on: 25.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "AdamicAdarIndex.h"

#include <cmath>

namespace NetworKit {

AdamicAdarIndex::AdamicAdarIndex(const Graph& G) : LinkPredictor(G) {
}

double AdamicAdarIndex::runImpl(node u, node v) {
  commonNeighborsIndex.setGraph(*G);
  std::vector<node> commonNeighbors = commonNeighborsIndex.getCommonNeighbors(u, v);
  double sum = 0;
  for (node w : commonNeighbors) {
    sum += 1.0 / std::log(G->degree(w));
  }
  return sum;
}


} // namespace NetworKit
