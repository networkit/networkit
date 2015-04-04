/*
 * LinkThresholder.h
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef LINKTHRESHOLDER_H_
#define LINKTHRESHOLDER_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Provides data points for the receiver operating characteristic of
 * a given set of predictions for graph edges.
 */
class LinkThresholder {
public:
  static std::vector<std::pair<node, node>> byScore(std::vector<LinkPredictor::node_dyad_score_pair> predictions, double minScore);

  static std::vector<std::pair<node, node>> byCount(std::vector<LinkPredictor::node_dyad_score_pair> predictions, count numLinks);
  
  static std::vector<std::pair<node, node>> byPercentage(std::vector<LinkPredictor::node_dyad_score_pair> predictions, double percentageLinks);
  
};

} // namespace NetworKit

#endif /* LINKTHRESHOLDER_H_std::vector<std::pair<node, node>> by */