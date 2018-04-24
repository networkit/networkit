/*
 * GroupDegree.cpp
 *
 *  Created on: 20.04.2018
 *      Author: Eugenio Angriman
 */

#include "GroupDegree.h"

namespace NetworKit {
GroupDegree::GroupDegree(const Graph &G, count k)
    : G(G), k(k), n(G.upperNodeIdBound()), queue(Aux::BucketPQ(n, -n + 1, 1)) {
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

  group.clear();
  group.reserve(k);
  inGroup.assign(n, false);

  reachable.assign(n, false);
  gain.assign(n, 0);
  hasComputedScore = false;
}

void GroupDegree::run() {
  init();

  G.forNodes([&](node u) {
    int64_t curNodeScore = G.degree(u);
    queue.insert(-curNodeScore, u);
    gain[u] = curNodeScore;
  });

  group.push_back(queue.extractMin().second);
  inGroup[group.back()] = true;
  while (group.size() < k) {
    updateQueue();
    group.push_back(queue.extractMin().second);
    inGroup[group.back()] = true;
  }

  hasRun = true;
}

void GroupDegree::updateQueue() {
  node lastAdded = group.back();
  INFO("Adding node ", lastAdded);
  reachable[lastAdded] = true;
  affected.assign(n, false);
  std::vector<node> neighbors = G.neighbors(lastAdded);
  if (G.isDirected()) {
    G.forInNeighborsOf(lastAdded, [&](node v) {
      affected[v] = true;
      INFO("Marked ", v, " as affected");
    });
  }
  //#pragma omp parallel for
  for (count i = 0; i < neighbors.size(); ++i) {
    node u = neighbors[i];
    if (!inGroup[u]) {
      INFO("Gain of neighbor ", u, " diminished from ", gain[u] + 1, " to ",
           gain[u]);
      reachable[u] = true;
      affected[u] = true;
      if (G.isDirected()) {
        G.forInNeighborsOf(u, [&](node v) {
          if (!affected[v] && !inGroup[v]) {
            //#pragma omp critical
            affected[v] = true;
            INFO("Marking ", v, " as affected");
          }
        });
      } else {
        G.forNeighborsOf(u, [&](node v) {
          if (!affected[v] && !inGroup[v]) {
            //#pragma omp critical
            affected[v] = true;
            INFO("Marking ", v, " as affected");
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
      if (G.isDirected()) {
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
        if (!G.isDirected()) {
          if (!groupNeighbor && inGroup[v]) {
            groupNeighbor = true;
            --newGain;
          }
        }
      });
      gain[i] = newGain;
      INFO("New gain of affected node ", i, " = ", newGain);
#pragma omp critical
      queue.changeKey(-newGain, i);
    }
  }
}
} // namespace NetworKit
