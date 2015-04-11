/*
 * AdamicAdarIndex.cpp
 *
 *  Created on: 25.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "AdamicAdarIndex.h"
#include "NeighborhoodUtility.h"

#include <cmath>

namespace NetworKit {

double AdamicAdarIndex::runImpl(node u, node v) {
  std::vector<node> commonNeighbors = NeighborhoodUtility::getCommonNeighbors(*G, u, v);
  double sum = 0;
  for (node w : commonNeighbors) {
    sum += 1.0 / std::log(G->degree(w));
  }
  return sum;
}

} // namespace NetworKit
