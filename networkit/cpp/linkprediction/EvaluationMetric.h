/*
 * EvaluationMetric.h
 *
 *  Created on: 17.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef EVALUATIONMETRIC_H_
#define EVALUATIONMETRIC_H_

#include "../graph/Graph.h"
#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 * Abstract base class for evaluation curves of link predictors.
 * The evualation curves are generated based on the node-dyad-score pairs calculated
 * by the link predictor and a testGraph used for performance evaluation.
 */
class EvaluationMetric {
private:

  /**
   * Implementation of the generatePoints-method.
   * Doesn't have to check if testGraph or predictions are existent.
   */
  virtual void generatePointsImpl() = 0;

  void calculateStatisticalMeasures();

  void addStatisticsForThresholdIndex(index thresholdIndex);

  void setPositivesAndNegatives();

  void setTrueAndFalsePositives(std::vector<index> thresholdIndices);
  
  void setTrueAndFalseNegatives(std::vector<index> thresholdIndices);

protected:
  std::pair<std::vector<double>, std::vector<double>> generatedPoints; //!< The generated points of the curve

  const Graph* testGraph; //!< Contains set of edges to test the given node-pairs and its scores against

  std::vector<LinkPredictor::node_dyad_score_pair> predictions; //!< Pairs of node-pairs and corresponding scores generated from the LinkPredictor to evaluate

  count numThresholds;

  count numPositives;

  count numNegatives;

  std::vector<count> truePositives;

  std::vector<count> falsePositives;

  std::vector<count> trueNegatives;

  std::vector<count> falseNegatives;

public:
  explicit EvaluationMetric();

  /**
   *
   * @param testGraph Graph containing test-set of edges to use for evaluation
   * @param predictions Dyad-score-pairs whose prediction quality has to be evaluated
   */
  explicit EvaluationMetric(const Graph& testGraph, 
      std::vector<LinkPredictor::node_dyad_score_pair> predictions = std::vector<LinkPredictor::node_dyad_score_pair>());

  virtual ~EvaluationMetric() = default;

  void setTestGraph(const Graph& newTestGraph);

  void setPredictions(std::vector<LinkPredictor::node_dyad_score_pair> newPredictions);

  /**
   * Generates the data-points for the EvaluationMetric.
   *
   * @param numThresholds Number of thresholds to use.
   * If the number exceeds the actual number of predictions all the predictions will be used.
   * If set to 0 the maximal number of thresholds will be used.
   */
  void generatePoints(count numThresholds = 1000);

  /**
   * Returns the previously generated data-points for the curve.
   * @return the previously generated data-points for the curve where the first vector
   * of the pair contains all x-values and the second vector the corresponding y-values
   */
  std::pair<std::vector<double>, std::vector<double>> getPoints() const;

  /**
   * Calculates the area under the curve for the previously generated data-points.
   * This is done through the trapezoid rule.
   * @return the area under the given curve
   */
  virtual double areaUnderCurve() const;

};

} // namespace NetworKit

#endif /* EVALUATIONMETRIC_H_ */