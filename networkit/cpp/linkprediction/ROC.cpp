/*
 * ROC.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "ROC.h"

namespace NetworKit {

void ROC::generatePoints() {
  std::vector<std::pair<std::pair<node, node>, double>> predictedEdges;
  generatedPoints.first.clear();
  generatedPoints.second.clear();
  count positives = 0;
  for (std::pair<std::pair<node, node>, double> p : dyadScorePairs) {
    if (testGraph.hasEdge(p.first.first, p.first.second)) {
      positives++;
    }
  }
  count negatives = dyadScorePairs.size() - positives;
  count truePositives = 0;
  count falsePositives = 0;
  for (index i = 0; i <= dyadScorePairs.size(); ++i) {
    truePositives = 0;
    falsePositives = 0;
    for (index j = 0; j < dyadScorePairs.size(); ++j) {
      if (j < i) { // predicted as an edge
        if (testGraph.hasEdge(dyadScorePairs[j].first.first, dyadScorePairs[j].first.second)) {
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