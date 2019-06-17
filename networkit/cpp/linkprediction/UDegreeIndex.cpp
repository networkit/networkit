/*
 * UDegreeIndex.cpp
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <networkit/linkprediction/UDegreeIndex.hpp>

namespace NetworKit {

double UDegreeIndex::runImpl(node u, node) {
  return G->degree(u);
}

} // namespace NetworKit
