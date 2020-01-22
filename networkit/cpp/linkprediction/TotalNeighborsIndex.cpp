/*
 * TotalNeighborsIndex.cpp
 *
 *  Created on: 05.04.2015
 *      Author: Kolja Esders
 */

#include <networkit/linkprediction/NeighborhoodUtility.hpp>
#include <networkit/linkprediction/TotalNeighborsIndex.hpp>

namespace NetworKit {

double TotalNeighborsIndex::runImpl(node u, node v) {
  return NeighborhoodUtility::getNeighborsUnion(*G, u, v).size();
}

} // namespace NetworKit
