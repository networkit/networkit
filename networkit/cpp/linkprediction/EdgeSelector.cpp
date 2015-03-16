/*
 * EdgeSelector.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "EdgeSelector.h"

#include "../auxiliary/Log.h"

namespace NetworKit {

EdgeSelector::EdgeSelector(const Graph& G, LinkPredictor* linkPredictor, count maxSize)
    : G(G), linkPredictor(linkPredictor), maxSize(maxSize) {
}

const std::vector<EdgeSelector::nodes_score_pair> EdgeSelector::getAll() const {
  if (!executed) {
    throw std::logic_error("Call calculateScores() first.");
  }
  return nodePairs;
}

std::vector<EdgeSelector::nodes_score_pair> EdgeSelector::getByCount(count cnt) const {
  if (maxSize > 0 && cnt > maxSize) {
    throw std::invalid_argument("Given count exceeds maximal number of elements to store set during construction.");
  } else if (!executed) {
    throw std::logic_error("Call calculateScores() first.");
  }
  std::vector<nodes_score_pair> limitedNodePairs(cnt);
  std::copy(nodePairs.begin(), nodePairs.begin() + cnt, limitedNodePairs.begin());
  return limitedNodePairs;
}

void EdgeSelector::calculateScores() {
  std::priority_queue<nodes_score_pair, std::vector<nodes_score_pair>, SecondGreater> pairQueue;
  G.forNodes([&](node u) {
    G.forNodes([&](node v) {
      if (u < v && !G.hasEdge(u, v)) {
        double score = linkPredictor->run(u, v);
        if (maxSize == 0 || pairQueue.size() < maxSize) {
          pairQueue.push(std::make_pair(std::make_pair(u, v), score));
        } else if (score > pairQueue.top().second) {
          pairQueue.pop();
          pairQueue.push(std::make_pair(std::make_pair(u, v), score));
        }
      }
    });
  });
  //INFO("SizeOf(pairQueue) = ", pairQueue.size());
  count numEntries = pairQueue.size();
  for (index i = 0; i < numEntries; ++i) {
    nodePairs.push_back(pairQueue.top());
    pairQueue.pop();
  }
  executed = true;
}

//std::vector<std::pair<node, node>> EdgeSelector::selectByThreshold(double threshold) {
  // Use it != topPredictions.lower_bound(10) in loop.
//}

} // namespace NetworKit