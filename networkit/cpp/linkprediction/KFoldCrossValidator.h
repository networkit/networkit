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
 * Averages the performance of a given LinkPredictor on a given graph.
 * This is done by randomly partitioning the graph into k subsamples.
 * From the k samples a single subsamples is used as test-data and the
 * remaining subsamples are used as training-data. This is done for every subsample.
 * The single performances will be measured by the given EvaluationMetric.
 */
class KFoldCrossValidator {
private:
  const Graph& G; //!< The graph used as testGraph

  LinkPredictor* linkPredictor; //!< Predictor whose performance to evaluate

  EvaluationMetric* evaluator; //!< Metric used to generate an evaluation score

public:
  /**
   *
   * @param G Graph to work on
   * @param linkPredictor Predictor whose performance should be measured
   * @param evaluator Evaluator which provides a metric for evaluating the link prediction results
   */
  explicit KFoldCrossValidator(const Graph& G, LinkPredictor* linkPredictor, EvaluationMetric* evaluator);
  
  /**
   * Calculates the average AUC of the given EvaluationMetric after @a k test-runs.
   * @param k Number of subsamples to split the given Graph G into
   * @return the average area under the given EvaluationMetric after @a k runs
   */
  double crossValidate(count k);
  
};

} // namespace NetworKit

#endif /* KFOLDCROSSVALIDATOR_H_ */