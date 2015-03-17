/*
 * EvaluationCurve.h
 *
 *  Created on: 17.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef EVALUATIONCURVE_H_
#define EVALUATIONCURVE_H_

#include "../graph/Graph.h"
#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 * Abstract base class for evaluation curves of link predictors.
 * The evualation curves are generated based on the node-dyad-score pairs calculated
 * by the link predictor and a testGraph used for performance evaluation.
 */
class EvaluationCurve {
protected:
  std::pair<std::vector<double>, std::vector<double>> generatedPoints; //!< The generated points of the curve

  const Graph& testGraph; //!< Contains set of edges to test the given node-pairs and its scores against

  std::vector<LinkPredictor::node_dyad_score_pair> dyadScorePairs; //!< Pairs of node-pairs and corresponding scores generated from the LinkPredictor to evaluate

public:
  /**
   *
   * @param testGraph Graph containing test-set of edges to use for evaluation
   * @param dyadScorePairs Dyad-score-pairs whose prediction quality has to be evaluated
   */
  EvaluationCurve(const Graph& testGraph, std::vector<LinkPredictor::node_dyad_score_pair> dyadScorePairs);

  virtual ~EvaluationCurve() = default;

  /**
   * Generates the data-points for the EvaluationCurve.
   */
  virtual void generatePoints() = 0;

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

#endif /* EVALUATIONCURVE_H_ */