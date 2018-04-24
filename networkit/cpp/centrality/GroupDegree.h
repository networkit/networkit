/*
 * GroupDegree.h
 *
 *  Created on: 20.04.2018
 *      Author: Eugenio Angriman
 */

#ifndef GROUPDEGREE_H_
#define GROUPDEGREE_H_

#include "../auxiliary/BucketPQ.h"
#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

class GroupDegree : public Algorithm {

public:
  GroupDegree(const Graph &G, count k = 1);
  void run() override;
  std::vector<node> groupMaxDegree();
  count getScore();

protected:
  Graph G;
  const count k;
  count n;
  std::vector<node> group;
  std::vector<int64_t> gain;
  std::vector<bool> reachable;
  std::vector<bool> affected;
  std::vector<bool> inGroup;
  Aux::BucketPQ queue;
  count groupScore;
  bool hasComputedScore;
  bool hasSortedGroup;

  void init();
  void updateQueue();
  void computeScore();
  void checkHasRun();
};

inline std::vector<node> GroupDegree::groupMaxDegree() {
  checkHasRun();
  if (!hasSortedGroup) {
    std::sort(group.begin(), group.end());
    hasSortedGroup = true;
  }
  return group;
}

inline count GroupDegree::getScore() {
  checkHasRun();
  if (!hasComputedScore) {
    computeScore();
  }
  return groupScore;
}

inline void GroupDegree::computeScore() {
  groupScore = 0;
  G.forNodes([&](node u) {
    if (reachable[u] && !inGroup[u]) {
      ++groupScore;
    }
  });
  hasComputedScore = true;
}

inline void GroupDegree::checkHasRun() {
  if (!hasRun) {
    throw std::runtime_error("Run method has not been called.");
  }
}
} // namespace NetworKit

#endif
