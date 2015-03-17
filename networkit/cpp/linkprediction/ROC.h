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
   * Generates a vector of points that belong to the Receiver Operating Characteristic of the given
   * dyad-score-pairs and the graph to check against.
   * @param testGraph The graph whose edges are used to test the given dyad-score-pairs against
   * @param data Dyad-score-pairs to test.
   * @return a pair of vectors where the first vector contains the x values and the second the y values
   */
  static std::pair<std::vector<double>, std::vector<double>> fromDyadScorePairs(const Graph &testGraph, std::vector<LinkPredictor::node_dyad_score_pair> data);
};

} // namespace NetworKit

#endif /* ROC_H_ */