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

std::vector<std::pair<double, double>> ROC::from(const Graph& testGraph, std::vector<LinkPredictor::node_dyad_score_pair> data) {
  std::vector<std::pair<std::pair<node, node>, double>> predictedEdges;
  std::vector<std::pair<double, double>> points;
  count positives = 0;
  for (std::pair<std::pair<node, node>, double> p : data) {
    if (testGraph.hasEdge(p.first.first, p.first.second)) {
      positives++;
    }
  }
  count negatives = data.size() - positives;
  INFO("Positives = ", positives);
  INFO("Negatives = ", negatives);
  INFO("Overall undirected edges = ", testGraph.numberOfEdges());
  count truePositives = 0;
  count falsePositives = 0;
  for (index i = 0; i <= data.size(); ++i) {
    INFO("Threshold: ", i, " / ", data.size());
    truePositives = 0;
    falsePositives = 0;
    for (index j = 0; j < data.size(); ++j) {
      if (j < i) { // predicted as an edge
        INFO("Predicted edge: (", data[j].first.first, ", ", data[j].first.second, ") with score ", data[j].second, ".");
        if (testGraph.hasEdge(data[j].first.first, data[j].first.second)) {
          truePositives++;
        } else {
          falsePositives++;
        }
      }
    }
    points.push_back(std::make_pair(1.0 * falsePositives / negatives, 1.0 * truePositives / positives));
    INFO("Generated point (", falsePositives," / ", negatives, ", ", truePositives, " / ", positives, ") [", points.back().first, ", ", points.back().second, "].\n");
  }
  return points;
}

} // namespace NetworKit