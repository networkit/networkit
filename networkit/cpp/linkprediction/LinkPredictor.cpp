/*
 * LinkPredictor.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <algorithm>

#include "LinkPredictor.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

LinkPredictor::LinkPredictor() : G(nullptr), validCache(false) {
}

LinkPredictor::LinkPredictor(const Graph& G) : G(&G), validCache(false) {
}

void LinkPredictor::setGraph(const Graph& newGraph) {
  G = &newGraph;
  validCache = false;
}

std::vector<LinkPredictor::node_dyad_score_pair> LinkPredictor::runOn(std::vector<std::pair<node, node>> nodePairs) {
  std::vector<node_dyad_score_pair> predictions;
  for (index i = 0; i < nodePairs.size(); ++i) {
    predictions.push_back(std::make_pair(nodePairs[i], run(nodePairs[i].first, nodePairs[i].second)));
  }
  std::sort(predictions.begin(), predictions.end(), ConcreteNodeDyadScoreComp);
  return predictions;
}

std::vector<LinkPredictor::node_dyad_score_pair> LinkPredictor::runOnParallel(std::vector<std::pair<node, node>> nodePairs) {
  std::vector<node_dyad_score_pair> predictions;
  #pragma omp parallel
  {
    std::vector<node_dyad_score_pair> predictionsPrivate;
    #pragma omp for nowait schedule(dynamic, 128)
    for (index i = 0; i < nodePairs.size(); ++i) {
      predictionsPrivate.push_back(std::make_pair(nodePairs[i], run(nodePairs[i].first, nodePairs[i].second)));
    }
    #pragma omp critical
    predictions.insert(predictions.end(), predictionsPrivate.begin(), predictionsPrivate.end());
  }
  std::sort(predictions.begin(), predictions.end(), ConcreteNodeDyadScoreComp);
  return predictions;
}

double LinkPredictor::run(node u, node v) {
  if (G == nullptr) {
    throw std::logic_error("Set a graph first.");
  } else if (!G->hasNode(u) || !G->hasNode(v)) {
    throw std::invalid_argument("Invalid node provided.");
  }
  return runImpl(u, v);
}

std::vector<LinkPredictor::node_dyad_score_pair> LinkPredictor::runAll() {
  std::vector<node> nodes = G->nodes();
  std::vector<std::pair<node, node>> nodePairs;
  for (index i = 0; i < nodes.size(); ++i) {
    for (index j = i + 1; j < nodes.size(); ++j) {
      if (!G->hasEdge(i, j))
        nodePairs.push_back(std::make_pair(i, j));
    }
  }
  return runOn(nodePairs);
}

} // namespace NetworKit