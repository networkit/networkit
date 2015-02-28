/*
 * EdgeSelector.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "EdgeSelector.h"

namespace NetworKit {

EdgeSelector::EdgeSelector(const Graph& G, const LinkPredictor linkPredictor) : G(G), linkPredictor(linkPredictor) {
}

std::vector<std::pair<node, node>> EdgeSelector::selectByLimit(count limit) {
  // TODO
  std::map<double, std::pair<node, node>> topPredictions;
  return new std::vector<std::pair<>>();
}

} // namespace NetworKit