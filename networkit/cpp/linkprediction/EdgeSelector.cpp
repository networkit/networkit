/*
 * EdgeSelector.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "EdgeSelector.h"

namespace NetworKit {

EdgeSelector::EdgeSelector(const Graph& G, LinkPredictor* linkPredictor) : G(G), linkPredictor(linkPredictor) {
}

std::vector<std::pair<node, node>> EdgeSelector::selectByLimit(count limit) const {
  std::map<double, std::pair<node, node>> topPredictions;
  G.forNodes([&](node u) {
    G.forNodes([&](node v) {
      if (u != v && !G.hasEdge(u, v)) {
        // TODO: Investigate parallel insert into multimap
        double score = linkPredictor->run(u, v);
        if (topPredictions.size() < limit) {
          topPredictions.insert(std::make_pair(score, std::make_pair(u, v)));
        } else if (score > (*topPredictions.begin()).first) {
          topPredictions.erase(topPredictions.begin());
          topPredictions.insert(std::make_pair(score, std::make_pair(u, v)));
        }
      }
    });
  });
  // Return edges with the highest scores
  std::vector<std::pair<node, node>> topEdges;
  for (std::map<double, std::pair<node, node>>::iterator it = topPredictions.begin(); it != topPredictions.end(); ++it) {
    topEdges.push_back((*it).second);
  }
  return topEdges;
}

//std::vector<std::pair<node, node>> EdgeSelector::selectByThreshold(double threshold) {
  // Use it != topPredictions.lower_bound(10) in loop.
//}

} // namespace NetworKit