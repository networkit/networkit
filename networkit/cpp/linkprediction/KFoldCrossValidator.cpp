/*
 * KFoldCrossValidator.cpp
 *
 *  Created on: 18.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "KFoldCrossValidator.h"
#include "RandomEdgePartitioner.h"

namespace NetworKit {

KFoldCrossValidator::KFoldCrossValidator(const Graph& G, LinkPredictor* linkPredictor, EvaluationMetric* evaluator)
    : G(G), linkPredictor(linkPredictor), evaluator(evaluator) {
}

double KFoldCrossValidator::crossValidate(count k) {
  if (k < 2) {
    throw std::invalid_argument("k < 2 => Unable to split set of edges.");
  }
  std::vector<Graph> edgeSets;
  count setSize = G.numberOfEdges() / k;
  Graph remainingGraph(G);
  // Split edge-set into k subsets
  for (count i = 0; i < k;  ++i) {
    RandomEdgePartitioner partitioner(remainingGraph);
    std::pair<Graph, Graph> splittedEdges = partitioner.partitionByCount(setSize);
    edgeSets.push_back(splittedEdges.second);
  }
  // perform evaluation k-times and average over result
  double sum = 0;
  for (count i = 0; i < k; ++i) {
    Graph testSet = edgeSets.at(i);
    Graph trainingSet(G);
    testSet.forEdges([&](node u, node v) {
      trainingSet.removeEdge(u, v);
    });
    linkPredictor->setGraph(trainingSet);
    evaluator->setTestGraph(testSet);
    evaluator->setPredictions(linkPredictor->runAll());
    evaluator->generatePoints();
    double auc = evaluator->areaUnderCurve();
    sum += auc;
  }
  return sum / k;
}



} // namespace NetworKit