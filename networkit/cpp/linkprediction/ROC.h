/*
 * ROC.h
 *
 *  Created on: 14.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef ROC_H_
#define ROC_H_

#include <memory>

#include "../graph/Graph.h"
#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Provides data points for the receiver operating characteristic of
 * a given set of predictions for graph edges.
 */
class ROC {
public:
  /**
   * Constructs a new LinkPredictor instance for the graph @a G.
   *
   * @param G The graph to use
   * @param predictions The edges which got predicted
   */
  //explicit ROC(const Graph& testGraph);
  
  /**
   * 
   *
   * @return a prediction-score indicating the likelihood of an
   * edge between the given nodes
   */
  static std::vector<std::pair<double, double>> from(const Graph &testGraph, std::vector<LinkPredictor::node_dyad_score_pair> data);
};

} // namespace NetworKit

#endif /* ROC_H_ */