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

  groupScore = 0;
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
  score.assign(n, 0);
}

void GroupDegree::run() {
  init();

  G.forNodes([&](node u) {
    int64_t curNodeScore = G.degree(u);
    queue.insert(-curNodeScore, u);
    score[u] = curNodeScore;
  });

  pushInQueue();
  while (group.size() < k) {
    updateQueue(group.back());
    pushInQueue();
  }

  hasRun = true;
}

void GroupDegree::updateQueue(node lastAdded) {
  reachable[lastAdded] = true;
  G.forNeighborsOf(lastAdded, [&](node u) { reachable[u] = true; });
  if (G.isDirected()) {

  } else {
    G.forNeighborsOf(lastAdded, [&](node u) {
      count newScore = 0;
      tmp.clear();
      tmp.reserve(n);
      G.forNeighborsOf(u, [&](node v) {
        if (!reachable[v]) {
          tmp.push_back(v);
          ++newScore;
        }
      });
#pragma omp parallel for
      for (count i = 0; i < tmp.size(); ++i) {
        count neighborNewScore = 0;
        node curNeighbor = tmp[i];
        G.forNeighborsOf(curNeighbor, [&](node v) {
          if (!reachable[v]) {
            ++neighborNewScore;
          }
        });
        score[curNeighbor] = neighborNewScore;
        queue.changeKey(-score[curNeighbor], curNeighbor);
      }
      score[u] = (newScore == 1) ? 0 : newScore;
      queue.changeKey(-score[u], u);
    });
  }
}

void GroupDegree::pushInQueue() {
  auto highestDegree = queue.extractMin();
  group.push_back(highestDegree.second);
  inGroup[highestDegree.second] = true;
  groupScore += -highestDegree.first;
}
} // namespace NetworKit
