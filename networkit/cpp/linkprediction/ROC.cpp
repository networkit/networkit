/*
 * ROC.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "ROC.h"

namespace NetworKit {

void ROC::generatePointsImpl() {
  std::vector<std::pair<std::pair<node, node>, double>> predictedEdges;
  generatedPoints.first.clear();
  generatedPoints.second.clear();
  count positives = 0;
  for (std::pair<std::pair<node, node>, double> p : predictions) {
    if (testGraph->hasEdge(p.first.first, p.first.second)) {
      positives++;
    }
  }
  count negatives = predictions.size() - positives;
  count truePositives = 0;
  count falsePositives = 0;
  for (index i = 0; i <= predictions.size(); ++i) {
    truePositives = 0;
    falsePositives = 0;
    for (index j = 0; j < predictions.size(); ++j) {
      if (j < i) { // predicted as an edge
        if (testGraph->hasEdge(predictions[j].first.first, predictions[j].first.second)) {
          truePositives++;
        } else {
          falsePositives++;
        }
      }
    }
    double falsePositiveRatio = 1.0 * falsePositives / negatives;
    if (generatedPoints.first.size() == 0 || generatedPoints.first.back() < falsePositiveRatio) {
      generatedPoints.first.push_back(falsePositiveRatio);
      generatedPoints.second.push_back(1.0 * truePositives / positives);
    }
  }
}

} // namespace NetworKit