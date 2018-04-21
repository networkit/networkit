/*
 * GroupDegree.cpp
 *
 *  Created on: 20.04.2018
 *      Author: Eugenio Angriman
 */

#include "GroupDegree.h"

namespace NetworKit {
GroupDegree::GroupDegree(const Graph &G, count k) : G(G), k(k) {
  if (k > G.upperNodeIdBound() || k <= 0) {
    throw std::runtime_error("k must be between 1 and n");
  }
}

void GroupDegree::init() {
  group.clear();
  group.reserve(k);

  n = G.upperNodeIdBound();
  inGroup.assign(n, false);
  score.assign(n, 0);
  G.removeSelfLoops();
}

void GroupDegree::run() {
  init();

  Aux::BucketPQ queue(n, -n + 1, 0);

  G.forNodes([&](node u) {
    queue.insert(-G.degree(u), u);
    score[u] = G.degree(u);
  });

  pushInQueue(queue);

  while (group.size() < k) {
    updateQueue(group.back(), queue);
    pushInQueue(queue);
  }

  hasRun = true;
}

void GroupDegree::pushInQueue(Aux::BucketPQ &queue) {
  auto highestDegree = queue.extractMin();
  group.push_back(highestDegree.second);
  inGroup[highestDegree.second] = true;

  groupScore += -highestDegree.first;
  for (node u : G.neighbors(highestDegree.second)) {
    if (inGroup[u]) {
      groupScore -= 1;
      return;
    }
  }
}

void GroupDegree::updateQueue(node lastAdded, Aux::BucketPQ &queue) {
  G.forNeighborsOf(lastAdded, [&](node u) {
    if (!inGroup[u]) {
      --score[u];
      queue.changeKey(-score[u], lastAdded);
    }
  });
}
} // namespace NetworKit
