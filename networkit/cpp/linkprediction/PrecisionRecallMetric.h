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
  std::pair<std::vector<double>, std::vector<double>> generatePoints() override;

public:
  using EvaluationMetric::EvaluationMetric;
  
};

} // namespace NetworKit

#endif /* PRECISIONRECALLMETRIC_H_ */