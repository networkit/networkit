/*
 * AdjustedRandIndex.cpp
 *
 *  Created on: 11.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "AdjustedRandIndex.h"

namespace NetworKit {

double AdjustedRandIndex::runImpl(node u, node v) {
  std::vector<node> uAdjacencyColumn(G->upperNodeIdBound());
  std::vector<node> vAdjacencyColumn(G->upperNodeIdBound());
  double a = 0, b = 0, c = 0, d = 0;
  G->forNeighborsOf(u, [&](node z) {
    uAdjacencyColumn[z] = 1;
  });
  G->forNeighborsOf(v, [&](node z) {
    vAdjacencyColumn[z] = 1;
  });
  for (index i = 0; i < G->upperNodeIdBound(); ++i) {
    a += uAdjacencyColumn[i] * vAdjacencyColumn[i];
    b += uAdjacencyColumn[i] * (1 - vAdjacencyColumn[i]);
    c += (1 - uAdjacencyColumn[i]) * vAdjacencyColumn[i];
    d += (1 - uAdjacencyColumn[i]) * (1 - vAdjacencyColumn[i]);
  }
  double ad = a*d;
  return (2*(ad - b*c)) / (a*b + a*c + 2*ad + b*b + b*d + c*c + c*d);
}

} // namespace NetworKit
