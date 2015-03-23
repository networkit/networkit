/*
 * LinkPredictor.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <algorithm>

#include "LinkPredictor.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

LinkPredictor::LinkPredictor() : G(nullptr), validCache(false) {
}

LinkPredictor::LinkPredictor(const Graph& G) : G(&G), validCache(false) {
}

void LinkPredictor::setGraph(const Graph& newGraph) {
  G = &newGraph;
  validCache = false;
}

std::vector<LinkPredictor::node_dyad_score_pair> LinkPredictor::runAllImpl(count limit) {
  std::priority_queue<node_dyad_score_pair, std::vector<node_dyad_score_pair>, NodeDyadScoreComp> pairQueue;
  std::vector<node_dyad_score_pair> result;
  G->forNodes([&](node u) {
    G->forNodes([&](node v) {
      if (u < v && !G->hasEdge(u, v)) {
        double score = run(u, v);
        if (limit == 0 || pairQueue.size() < limit) {
          //INFO("Inserting: ((", u, ", ", v, "), ", score, ")");
          pairQueue.push(std::make_pair(std::make_pair(u, v), score));
        } else if (score > pairQueue.top().second) {
          //INFO("Replacing ((", pairQueue.top().first.first, ", ", pairQueue.top().first.second, "), ", pairQueue.top().second, ") with ((", u, ", ", v, "), ", score, ")");
          pairQueue.pop();
          pairQueue.push(std::make_pair(std::make_pair(u, v), score));
        }
        //INFO("Top element: ((", pairQueue.top().first.first, ", ", pairQueue.top().first.second, "), ", pairQueue.top().second, ")");
      }
    });
  });
  count numEntries = pairQueue.size();
  for (index i = 0; i < numEntries; ++i) {
    result.push_back(pairQueue.top());
    pairQueue.pop();
  }
  std::reverse(result.begin(), result.end());
  return result;
}

double LinkPredictor::run(node u, node v) {
  if (G == nullptr) {
    throw std::logic_error("Set a graph first.");
  } else if (!G->hasNode(u) || !G->hasNode(v)) {
    throw std::invalid_argument("Invalid node provided.");
  }
  return runImpl(u, v);
}

std::vector<LinkPredictor::node_dyad_score_pair> LinkPredictor::runAll(count limit) {
  if (G == nullptr) {
    throw std::logic_error("Set a graph first.");
  }
  return runAllImpl(limit);
}

} // namespace NetworKit