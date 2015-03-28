/*
 * KatzIndex.cpp
 *
 *  Created on: 30.01.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <list>

#include "KatzIndex.h"

namespace NetworKit {

KatzIndex::KatzIndex(count maxPathLength, double dampingValue)
    : maxPathLength(maxPathLength), dampingValue(dampingValue) {
}

KatzIndex::KatzIndex(const Graph& G, count maxPathLength, double dampingValue)
    : LinkPredictor(G), maxPathLength(maxPathLength), dampingValue(dampingValue) {
}

double KatzIndex::getScore(node u, node v) const {
  node endNode = lastStartNode == u ? v : u;
  if (lastScores.find(endNode) == lastScores.end()) {
    return 0;
  }
  return lastScores.at(endNode);
}

double KatzIndex::runImpl(node u, node v) {
  if (validCache && (lastStartNode == u || lastStartNode == v)) {
    return getScore(u, v);
  }
  lastScores.clear();
  validCache = false;
  std::list<node> toProcess;
  // Start at the node with less neighbors to potentially increase performance
  lastStartNode = G->degree(u) > G->degree(v) ? v : u;
  toProcess.push_front(lastStartNode);
  for (unsigned int pathLength = 1; pathLength <= maxPathLength; ++pathLength) {
    std::unordered_map<node, unsigned int> hits;
    for (std::list<node>::const_iterator it = toProcess.begin(); it != toProcess.end(); ++it) {
      const node current = *it;
      // TODO: Parallelize this part.
      G->forNeighborsOf(current, [&](node neighbor) {
        hits[neighbor] += 1;
      });
    }
    // Add found nodes to the todo-list for the next round and update scores
    double dampingFactor = pow(dampingValue, pathLength);
    toProcess.clear();
    for (auto kv : hits) {
      lastScores[kv.first] += dampingFactor * kv.second;
      toProcess.push_front(kv.first);
    }
  }
  validCache = true;
  return getScore(u, v);
}

std::vector<LinkPredictor::node_dyad_score_pair> KatzIndex::runOnParallel(std::vector<std::pair<node, node>> nodePairs) {
  std::vector<node_dyad_score_pair> predictions;
  #pragma omp parallel
  {
    // Create local KatzIndex
    KatzIndex katz(*G, maxPathLength, dampingValue);
    std::vector<node_dyad_score_pair> predictionsPrivate;
    #pragma omp for nowait
    for (index i = 0; i < nodePairs.size(); ++i) {
      predictionsPrivate.push_back(std::make_pair(nodePairs[i], katz.run(nodePairs[i].first, nodePairs[i].second)));
    }
    #pragma omp critical
    predictions.insert(predictions.end(), predictionsPrivate.begin(), predictionsPrivate.end());
  }
  std::sort(predictions.begin(), predictions.end(), ConcreteNodeDyadScoreComp);
  return predictions;
}

} // namespace NetworKit