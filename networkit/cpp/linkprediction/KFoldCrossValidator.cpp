/*
 * KFoldCrossValidator.cpp
 *
 *  Created on: 18.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "KFoldCrossValidator.h"
#include "MissingLinksFinder.h"

namespace NetworKit {

KFoldCrossValidator::KFoldCrossValidator(const Graph& G, LinkPredictor* linkPredictor, EvaluationMetric* evaluator)
    : G(G), linkPredictor(linkPredictor), evaluator(evaluator) {
}

double KFoldCrossValidator::crossValidate(count k) {
  if (k < 2) {
    throw std::invalid_argument("k has to be at least 2.");
  }
  evaluator->setTestGraph(G);
  std::vector<std::vector<std::pair<node, node>>> linkPartitions(k, std::vector<std::pair<node, node>>());
  count numberEdges = G.numberOfEdges();
  Graph copyOfG(G);
  // Partition the graph into k sets of edges
  for (index i = 0; i < numberEdges; ++i) {
    std::pair<node, node> e = copyOfG.randomEdge();
    copyOfG.removeEdge(e.first, e.second);
    linkPartitions[i % k].push_back(e);
  }
  // Remove edges from original graph for every edge-set to get the trainingGraph
  double sum = 0;
  for (index i = 0; i < k; ++i) {
    Graph trainingGraph(G);
    for (index j = 0; j < linkPartitions[i].size(); ++j) {
      trainingGraph.removeEdge(linkPartitions[i][j].first, linkPartitions[i][j].second);
    }
    linkPredictor->setGraph(trainingGraph);
    // Evaluate on all 2-hop missing links
    std::vector<std::pair<node, node>> twoHopsMissing = MissingLinksFinder(trainingGraph).findAtDistance(2);
    std::pair<std::vector<double>, std::vector<double>> ps = evaluator->getCurve(linkPredictor->runOnParallel(twoHopsMissing));
    sum += evaluator->getAreaUnderCurve();
  }
  return sum / k;
}



} // namespace NetworKit