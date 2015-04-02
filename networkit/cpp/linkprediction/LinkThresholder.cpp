/*
 * LinkThresholder.cpp
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "LinkThresholder.h"

namespace NetworKit {

std::vector<std::pair<node, node>> LinkThresholder::byScore(std::vector<LinkPredictor::node_dyad_score_pair> predictions, double minScore) {
  std::list<LinkPredictor::node_dyad_score_pair> predictionsList;
  std::copy(predictions.begin(), predictions.end(), std::back_inserter(predictionsList));
  predictionsList.erase(std::remove_if(predictionsList.begin(), predictionsList.end(),
      [&](const std::pair<std::pair<node, node>, double>& p) { return p.second < minScore; }), predictionsList.end());
  std::vector<std::pair<node, node>> selectedLinks;
  selectedLinks.reserve(predictionsList.size());
  std::transform(predictionsList.begin(), predictionsList.end(), std::back_inserter(selectedLinks),
      [](const LinkPredictor::node_dyad_score_pair& p){ return p.first; });
  std::sort(selectedLinks.begin(), selectedLinks.end());
  return selectedLinks;
}

std::vector<std::pair<node, node>> LinkThresholder::byCount(std::vector<LinkPredictor::node_dyad_score_pair> predictions, count numLinks) {
  if (numLinks > predictions.size()) {
    throw std::invalid_argument("numLinks > predictions.size().");
  }
  LinkPredictor::sortByScore(predictions);
  std::vector<LinkPredictor::node_dyad_score_pair> selectedPredictions(predictions.begin(), predictions.begin() + numLinks);
  std::vector<std::pair<node, node>> selectedLinks;
  selectedLinks.reserve(selectedPredictions.size());
  std::transform(selectedPredictions.begin(), selectedPredictions.end(), std::back_inserter(selectedLinks),
      [](const LinkPredictor::node_dyad_score_pair& p){ return p.first; });
  std::sort(selectedLinks.begin(), selectedLinks.end());
  return selectedLinks;
}

std::vector<std::pair<node, node>> LinkThresholder::byPercentage(std::vector<LinkPredictor::node_dyad_score_pair> predictions, double percentageLinks) {
  if (percentageLinks < 0 || percentageLinks > 1) {
    throw std::invalid_argument("Given percentage is not in [0, 1].");
  }
  return byCount(predictions, percentageLinks * predictions.size());
}

} // namespace NetworKit