/*
 * TotalNeighborsIndex.cpp
 *
 *  Created on: 05.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "TotalNeighborsIndex.h"
#include "NeighborhoodUtility.h"

namespace NetworKit {

double TotalNeighborsIndex::runImpl(node u, node v) {
  return NeighborhoodUtility::getNeighborsUnion(*G, u, v).size();
}

} // namespace NetworKit
