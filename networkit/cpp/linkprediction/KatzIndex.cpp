/*
 * KatzIndex.cpp
 *
 *  Created on: 30.01.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <list>

#include "KatzIndex.h"

namespace NetworKit {

KatzIndex::KatzIndex(const Graph& G, count maxPathLength, double dampingValue)
    : LinkPredictor(G), maxPathLength(maxPathLength), dampingValue(dampingValue) {
  // Invalidate start node in case the algorithm starts at node 0 the first time
  lastStartNode = -1;
}

double KatzIndex::run(node u, node v) {
  if (lastStartNode == u || lastStartNode == v) {
    return getScore(u, v);
  }
  std::list<node> toProcess;
  // Start at the node with less neighbors to potentially increase performance
  lastStartNode = G.degree(u) > G.degree(v) ? v : u;
  toProcess.push_front(lastStartNode);
  for (unsigned int pathLength = 1; pathLength <= maxPathLength; ++pathLength) {
    std::unordered_map<node, unsigned int> hits;
    for (std::list<node>::const_iterator it = toProcess.begin(); it != toProcess.end(); ++it) {
      const node current = *it;
      // TODO: Parallelize this part.
      G.forNeighborsOf(current, [&](node neighbor) {
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
  return getScore(u, v);
}

double KatzIndex::getScore(node u, node v) const {
  node endNode = lastStartNode == u ? v : u;
  if (lastScores.find(endNode) == lastScores.end()) {
    return 0;
  }
  return lastScores.at(endNode);
}

} // namespace NetworKit