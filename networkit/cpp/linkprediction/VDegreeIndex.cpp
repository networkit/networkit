/*
 * VDegreeIndex.cpp
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "VDegreeIndex.h"

namespace NetworKit {

VDegreeIndex::VDegreeIndex() {
}

VDegreeIndex::VDegreeIndex(const Graph& G) : LinkPredictor(G) {
}

double VDegreeIndex::runImpl(node u, node v) {
  return G->degree(v);
}

} // namespace NetworKit