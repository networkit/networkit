/*
 * UDegreeIndex.cpp
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "UDegreeIndex.h"

namespace NetworKit {

double UDegreeIndex::runImpl(node u, node v) {
  return G->degree(u);
}

} // namespace NetworKit