/*
 * PrecisionRecallMetric.h
 *
 *  Created on: 21.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef PRECISIONRECALLMETRIC_H_
#define PRECISIONRECALLMETRIC_H_

#include "../graph/Graph.h"
#include "EvaluationMetric.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 */
class PrecisionRecallMetric : public EvaluationMetric {
private:
  /**
   * Generates the points for the Precision-Recall curve with respect to the given predictions.
   * The curve assigns every recall-value a corresponding precision as the y-value.
   * In case of a tie regarding multiple y-values for a x-value the smallest (= latest) y-value will be used.
   * @return a pair of vectors where the first vector contains all recall-values and the second vector
   * the correspoding precision-values
   */
  std::pair<std::vector<double>, std::vector<double>> generatePoints() override;

public:
  using EvaluationMetric::EvaluationMetric;
  
};

} // namespace NetworKit

#endif /* PRECISIONRECALLMETRIC_H_ */