/*
 * ROC.h
 *
 *  Created on: 14.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef ROC_H_
#define ROC_H_

#include "../graph/Graph.h"
#include "EvaluationCurve.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Provides data points for the receiver operating characteristic of
 * a given set of predictions for graph edges.
 */
class ROC : public EvaluationCurve {
public:
  using EvaluationCurve::EvaluationCurve;

  void generatePoints() override;
};

} // namespace NetworKit

#endif /* ROC_H_ */