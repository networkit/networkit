/*
 * KFoldCrossValidator.h
 *
 *  Created on: 18.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef KFOLDCROSSVALIDATOR_H_
#define KFOLDCROSSVALIDATOR_H_

#include "../graph/Graph.h"
#include "LinkPredictor.h"
#include "EvaluationMetric.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Evaluates the performance of a given LinkPredictor on a given graph
 * by randomly partitioning the graph into k subsamples. Of the k samples
 * a single subsamples is used as test-data and the remaining subsamples are
 * used as training-data.
 * The performance will be measured by the given EvaluationMetric.
 */
class KFoldCrossValidator {
private:
  const Graph& G;

  LinkPredictor* linkPredictor;

  EvaluationMetric* evaluator;

  Graph mergeEdges(std::set<Graph> edgeSets) const;

public:
  /**
   *
   * @param G Graph to work on
   * @param linkPredictor Predictor whose performance should be measured
   * @param evaluator Evaluator which provides a metric for evaluating the link prediction results
   */
  KFoldCrossValidator(const Graph& G, LinkPredictor* linkPredictor, EvaluationMetric* evaluator);
  
  /**
   * Calculates the average AUC of the given EvaluationMetric after k test-runs.
   *
   * @param k Number of subsamples to split the given Graph G into
   * @return the average area under the given EvaluationMetric after k runs
   */
  double crossValidate(count k);

};

} // namespace NetworKit

#endif /* KFOLDCROSSVALIDATOR_H_ */