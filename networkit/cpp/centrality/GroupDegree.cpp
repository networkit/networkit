/*
 * GroupDegree.cpp
 *
 *  Created on: 20.04.2018
 *      Author: Eugenio Angriman
 */

#include "GroupDegree.h"

namespace NetworKit {
GroupDegree::GroupDegree(const Graph &G, count k, bool countGroupNodes)
    : G(G), k(k), countGroupNodes(countGroupNodes), n(G.upperNodeIdBound()),
      queue(Aux::BucketPQ(n, -n + 1, countGroupNodes ? 0 : 1)) {
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

  if (hasRun) {
    n = G.upperNodeIdBound();
    while (queue.size() > 0) {
      queue.extractMin();
    }

    hasRun = false;
  }

  hasSortedGroup = false;
  hasComputedScore = false;
  group.clear();
  group.reserve(k);
  inGroup.assign(n, false);
  reachable.assign(n, false);
  gain.assign(n, 0);
}

void GroupDegree::run() {
  init();

  G.forNodes([&](node u) {
    int64_t curNodeScore = G.degreeOut(u);
    if (countGroupNodes) {
      ++curNodeScore;
    }
    queue.insert(-curNodeScore, u);
    gain[u] = curNodeScore;
  });

  updateGroup();
  while (group.size() < k) {
    updateQueue();
    updateGroup();
  }

  hasRun = true;
}

void GroupDegree::updateGroup() {
  group.push_back(queue.extractMin().second);
  inGroup[group.back()] = true;
  reachable[group.back()] = true;
}

void GroupDegree::updateQueue() {
  node lastAdded = group.back();
  affected.assign(n, false);
  std::vector<node> neighbors = G.neighbors(lastAdded);

#pragma omp parallel for
  for (count i = 0; i < neighbors.size(); ++i) {
    node u = neighbors[i];
    if (!inGroup[u] && !reachable[u]) {
      affected[u] = true;
      reachable[u] = true;
      if (G.isDirected()) {
        G.forInNeighborsOf(u, [&](node v) {
          if (!affected[v] && !inGroup[v]) {
#pragma omp critical
            affected[v] = true;
          }
        });
      } else {
        G.forNeighborsOf(u, [&](node v) {
          if (!affected[v] && !inGroup[v]) {
#pragma omp critical
            affected[v] = true;
          }
        });
      }
    }
  }

#pragma omp parallel for
  for (count i = 0; i < n; ++i) {
    if (affected[i]) {
      int64_t newGain = 0;
      bool groupNeighbor = false;
      if (!countGroupNodes && G.isDirected()) {
        G.forInNeighborsOf(i, [&](node v) {
          if (!groupNeighbor && inGroup[v]) {
            newGain = -1;
            groupNeighbor = true;
          }
        });
      }
      G.forNeighborsOf(i, [&](node v) {
        if (!reachable[v]) {
          ++newGain;
        }
        if (!countGroupNodes && !G.isDirected()) {
          if (!groupNeighbor && inGroup[v]) {
            groupNeighbor = true;
            --newGain;
          }
        }
      });
      gain[i] = newGain;
#pragma omp critical
      queue.changeKey(-newGain, i);
    }
  }
}
} // namespace NetworKit
