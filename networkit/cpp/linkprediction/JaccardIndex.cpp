/*
 * JaccardIndex.cpp
 *
 *  Created on: 23.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "JaccardIndex.h"
#include "NeighborhoodUtility.h"

namespace NetworKit {

double JaccardIndex::runImpl(node u, node v) {
  count unionSize = NeighborhoodUtility::getNeighborsUnion(*G, u, v).size();
  if (unionSize == 0) {
    return 0;
  }
  return 1.0 * NeighborhoodUtility::getCommonNeighbors(*G, u, v).size() / unionSize;
}

} // namespace NetworKit
