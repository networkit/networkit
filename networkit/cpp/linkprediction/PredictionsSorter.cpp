/*
 * PredictionsSorter.cpp
 *
 *  Created on: 26.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "PredictionsSorter.h"
#include "../auxiliary/Parallel.h"

namespace NetworKit {

void PredictionsSorter::sortByScore(std::vector<LinkPredictor::prediction>& predictions) {
  Aux::Parallel::sort(predictions.begin(), predictions.end(), ConcreteScoreComp);
}

void PredictionsSorter::sortByNodePair(std::vector<LinkPredictor::prediction>& predictions) {
  Aux::Parallel::sort(predictions.begin(), predictions.end(), ConcreteNodePairComp);
}

} // namespace NetworKit