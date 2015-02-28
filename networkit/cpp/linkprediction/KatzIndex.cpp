/*
 * KatzIndex.cpp
 *
 *  Created on: 30.01.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <unordered_map>
#include <list>

#include "KatzIndex.h"

namespace NetworKit {

KatzIndex::KatzIndex(const Graph& G, unsigned int maxPathLength, double dampingValue)
    : LinkPredictor(G), maxPathLength(maxPathLength), dampingValue(dampingValue) {
}

double KatzIndex::run(node u, node v) {
  if (lastStartNode == u || lastStartNode == v) {
    return lastScores.at(lastStartNode == u ? v : u);
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
  return lastScores.at(lastStartNode == u ? v : u);
}

} // namespace NetworKit