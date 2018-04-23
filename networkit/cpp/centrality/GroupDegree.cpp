/*
 * GroupDegree.cpp
 *
 *  Created on: 20.04.2018
 *      Author: Eugenio Angriman
 */

#include "GroupDegree.h"

namespace NetworKit {
GroupDegree::GroupDegree(const Graph &G, count k)
    : G(G), k(k), n(G.upperNodeIdBound()), queue(Aux::BucketPQ(n, -n + 1, 0)) {
  if (k > G.upperNodeIdBound() || k <= 0) {
    throw std::runtime_error("k must be between 1 and n");
  }
  if (G.numberOfSelfLoops() > 0) {
    throw std::runtime_error(
        "Group degree does not support graphs with self loops. Call "
        "removeSelfLoops() to remove self loops from the graph.");
  }
}

void GroupDegree::init() {

  hasComputedScore = false;
  if (hasRun) {
    n = G.upperNodeIdBound();
    while (queue.size() > 0) {
      queue.extractMin();
    }
    hasRun = false;
  }

  group.clear();
  group.reserve(k);

  reachable.assign(n, false);
  gain.assign(n, 0);
}

void GroupDegree::run() {
  init();

  G.forNodes([&](node u) {
    int64_t curNodeScore = G.degree(u);
    if (curNodeScore > 0) {
      --curNodeScore;
    }
    queue.insert(-curNodeScore, u);
    gain[u] = curNodeScore;
  });

  group.push_back(queue.extractMin().second);
  while (group.size() < k) {
    updateQueue();
    group.push_back(queue.extractMin().second);
  }

  hasRun = true;
}

void GroupDegree::updateQueue() {
  node lastAdded = group.back();
  reachable[lastAdded] = true;
  std::vector<node> neighbors = G.neighbors(lastAdded);
  affected.assign(n, false);
#pragma omp parallel for
  for (count i = 0; i < neighbors.size(); ++i) {
    node u = neighbors[i];
    reachable[u] = true;
    affected[u] = true;
    G.forNeighborsOf(u, [&](node v) {
#pragma omp critical
      affected[v] = true;
    });
  }

#pragma omp parallel for
  for (count i = 0; i < n; ++i) {
    if (affected[i]) {
      int64_t newGain = -1;
      G.forNeighborsOf(i, [&](node v) {
        if (!reachable[v]) {
          ++newGain;
        }
      });
      gain[i] = newGain;
#pragma omp critical
      queue.changeKey(-newGain, i);
    }
  }
}
} // namespace NetworKit
