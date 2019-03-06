/*
 * PredictionsSorter.cpp
 *
 *  Created on: 26.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "../../include/networkit/linkprediction/PredictionsSorter.hpp"
#include "../../include/networkit/auxiliary/Parallel.hpp"

namespace NetworKit {

void PredictionsSorter::sortByScore(std::vector<LinkPredictor::prediction>& predictions) {
  Aux::Parallel::sort(predictions.begin(), predictions.end(), ConcreteScoreComp);
}

void PredictionsSorter::sortByNodePair(std::vector<LinkPredictor::prediction>& predictions) {
  Aux::Parallel::sort(predictions.begin(), predictions.end(), ConcreteNodePairComp);
}

PredictionsSorter::ScoreComp PredictionsSorter::ConcreteScoreComp{};
PredictionsSorter::NodePairComp PredictionsSorter::ConcreteNodePairComp{};

} // namespace NetworKit
