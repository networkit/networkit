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
#include "EvaluationCurve.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 */
class KFoldCrossValidator {
private:
  const Graph& G;

  LinkPredictor* linkPredictor;

  EvaluationCurve* evaluator;

  Graph mergeEdges(std::set<Graph> edgeSets) const;

public:
  /**
   *
   * @param G 
   * @param linkPredictor 
   * @param evaluator 
   */
  KFoldCrossValidator(const Graph& G, LinkPredictor* linkPredictor, EvaluationCurve* evaluator);
  
  double crossValidate(count k);

};

} // namespace NetworKit

#endif /* KFOLDCROSSVALIDATOR_H_ */