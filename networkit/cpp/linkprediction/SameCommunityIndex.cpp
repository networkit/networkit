/*
 * SameCommunityIndex.cpp
 *
 *  Created on: 07.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "../../include/networkit/linkprediction/SameCommunityIndex.hpp"
#include "../../include/networkit/community/PLM.hpp"

namespace NetworKit {

SameCommunityIndex::SameCommunityIndex() : LinkPredictor(), communities() {
}

SameCommunityIndex::SameCommunityIndex(const Graph& graph) : LinkPredictor(graph) {
  PLM communityDetection(graph);
  communityDetection.run();
  communities = communityDetection.getPartition();
}

void SameCommunityIndex::setGraph(const Graph& newGraph) {
  LinkPredictor::setGraph(newGraph);
  PLM communityDetection(newGraph);
  communityDetection.run();
  communities = communityDetection.getPartition();
}

double SameCommunityIndex::runImpl(node u, node v) {
  return communities[u] == communities[v] ? 1 : 0;
}

} // namespace NetworKit
