/*
 * PrecisionRecallMetric.cpp
 *
 *  Created on: 21.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "PrecisionRecallMetric.h"

namespace NetworKit {

void PrecisionRecallMetric::generatePointsImpl() {
  std::vector<std::pair<std::pair<node, node>, double>> predictedEdges;
  generatedPoints.first.clear();
  generatedPoints.second.clear();
  count truePositives = 0;
  count falsePositives = 0;
  count falseNegatives = 0;
  for (index i = 0; i <= predictions.size(); ++i) {
    truePositives = 0;
    falsePositives = 0;
    falseNegatives = 0;
    for (index j = 0; j < predictions.size(); ++j) {
      bool actualEdge = testGraph->hasEdge(predictions[j].first.first, predictions[j].first.second);
      bool predictedEdge = j < i;
      if (predictedEdge && actualEdge) {
        truePositives++;
      } else if (predictedEdge && !actualEdge) {
        falsePositives++;
      } else if (!predictedEdge && actualEdge) {
        falseNegatives++;
      }
    }
    // Make sure not to divide by 0
    double recall = truePositives == 0 ? 0 : 1.0 * truePositives / (truePositives + falseNegatives);
    double precision = truePositives == 0 ? 0 : 1.0 * truePositives / (truePositives + falsePositives);
    if (generatedPoints.first.size() == 0 || generatedPoints.first.back() < recall) {
      generatedPoints.first.push_back(recall);
      generatedPoints.second.push_back(precision);
    }
  }
}

} // namespace NetworKit