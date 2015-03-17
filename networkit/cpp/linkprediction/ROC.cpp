/*
 * ROC.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "../auxiliary/Log.h"

#include "ROC.h"

namespace NetworKit {

/*ROC::ROC(const Graph& testGraph) : testGraph(testGraph) {
}*/

std::pair<std::vector<double>, std::vector<double>> ROC::fromDyadScorePairs(const Graph& testGraph, std::vector<LinkPredictor::node_dyad_score_pair> data) {
  std::vector<std::pair<std::pair<node, node>, double>> predictedEdges;
  std::pair<std::vector<double>, std::vector<double>> points;
  count positives = 0;
  for (std::pair<std::pair<node, node>, double> p : data) {
    if (testGraph.hasEdge(p.first.first, p.first.second)) {
      positives++;
    }
  }
  count negatives = data.size() - positives;
  count truePositives = 0;
  count falsePositives = 0;
  double lastFalsePositiveRatio = -1;
  for (index i = 0; i <= data.size(); ++i) {
    truePositives = 0;
    falsePositives = 0;
    for (index j = 0; j < data.size(); ++j) {
      if (j < i) { // predicted as an edge
        if (testGraph.hasEdge(data[j].first.first, data[j].first.second)) {
          truePositives++;
        } else {
          falsePositives++;
        }
      }
    }
    double falsePositiveRatio = 1.0 * falsePositives / negatives;
    if (lastFalsePositiveRatio == -1 || lastFalsePositiveRatio < falsePositiveRatio) {
      //points.push_back(std::make_pair(falsePositiveRatio, 1.0 * truePositives / positives));
      points.first.push_back(falsePositiveRatio);
      points.second.push_back(1.0 * truePositives / positives);
      lastFalsePositiveRatio = falsePositiveRatio;
    }
  }
  return points;
}

} // namespace NetworKit