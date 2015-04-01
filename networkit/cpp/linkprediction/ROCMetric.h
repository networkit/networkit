/*
 * ROCMetric.h
 *
 *  Created on: 14.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef ROCMETRIC_H_
#define ROCMETRIC_H_

#include "../graph/Graph.h"
#include "EvaluationMetric.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Provides data points for the receiver operating characteristic of
 * a given set of predictions for graph edges.
 */
class ROCMetric : public EvaluationMetric {
private:
  std::pair<std::vector<double>, std::vector<double>> generatePoints() override;

public:
  using EvaluationMetric::EvaluationMetric;
  
};

} // namespace NetworKit

#endif /* ROCMETRIC_H_ */