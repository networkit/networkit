/*
 * NeighborsMeasureIndex.cpp
 *
 *  Created on: 05.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "NeighborsMeasureIndex.h"

namespace NetworKit {

double NeighborsMeasureIndex::runImpl(node u, node v) {
  double neighborConnections = 0;
  G->forNeighborsOf(u, [&](node uNeighbor) {
    G->forNeighborsOf(v, [&](node vNeighbor) {
      // Don't count self-loops
      if (uNeighbor == vNeighbor || G->hasEdge(uNeighbor, vNeighbor)) {
        ++neighborConnections;
      }
    });
  }); 
  return neighborConnections;
}

} // namespace NetworKit
