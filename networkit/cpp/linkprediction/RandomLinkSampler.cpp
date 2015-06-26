/*
 * RandomLinkSampler.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "RandomLinkSampler.h"

namespace NetworKit {

namespace RandomLinkSampler {

Graph byPercentage(const Graph& G, double trainPercentage) {
  if (trainPercentage < 0 || trainPercentage > 1) {
    throw std::invalid_argument("Given percentage is not in [0, 1].");
  }
  return byCount(G, trainPercentage * G.numberOfEdges());
}

Graph byCount(const Graph& G, count numTrainLinks) {
  if (numTrainLinks > G.numberOfEdges()) {
    throw std::invalid_argument("numTrainLinks > G.numberOfEdges().");
  }
  Graph trainingGraph(G);
  for (count i = 0; i < G.numberOfEdges() - numTrainLinks; ++i) {
    std::pair<node, node> edgeToRemove = trainingGraph.randomEdge();
    trainingGraph.removeEdge(edgeToRemove.first, edgeToRemove.second);
  }
  return trainingGraph;
}

} // namespace RandomLinkSampler

} // namespace NetworKit