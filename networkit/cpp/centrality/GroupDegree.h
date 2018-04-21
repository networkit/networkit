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
  std::vector<int64_t> score;
  std::vector<bool> inGroup;
  count groupScore;

  void init();
  void pushInQueue(Aux::BucketPQ &queue);
  void updateQueue(node lastAdded, Aux::BucketPQ &queue);

  void checkHasRun();
};

inline std::vector<node> GroupDegree::groupMaxDegree() {
  checkHasRun();
  return group;
}

inline void GroupDegree::checkHasRun() {
  if (!hasRun) {
    throw std::runtime_error("Run method has not been called.");
  }
}

inline count GroupDegree::getScore() {
  checkHasRun();
  return groupScore;
}
} // namespace NetworKit

#endif
