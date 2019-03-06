/*
 * CommonNeighborsIndex.cpp
 *
 *  Created on: 06.12.2014
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "../../include/networkit/linkprediction/CommonNeighborsIndex.hpp"
#include "../../include/networkit/linkprediction/NeighborhoodUtility.hpp"

namespace NetworKit {

double CommonNeighborsIndex::runImpl(node u, node v) {
  return NeighborhoodUtility::getCommonNeighbors(*G, u, v).size();
}

} // namespace NetworKit