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

/**
 * @ingroup centrality
 */
class GroupDegree : public Algorithm {

public:
  /**
   * Finds the group with the highest group degree centrality according to the
   * definition proposed in 'The centrality of groups and classes' by Everett et
   * al. (The Journal of mathematical sociology, 1999). This is a submodular but
   * non monotone function so the algorithm can find a solution that is at least
   * 1/2 of the optimum. Worst-case running time is quadratic, but usually
   * faster in real-world networks.
   * The 'countGroupNodes' option also count the nodes inside the group in the
   * score, this make the group degree monotone and submodular and the algorithm
   * is guaranteed to return a (1 - 1/e)-approximation of the optimal solution.
   *
   * @param G A graph.
   * @param k Size of the group of nodes
   * @param countGroupNodes if nodes inside the group should be counted in the
   * centrality score.
   */
  GroupDegree(const Graph &G, count k = 1, bool countGroupNodes = false);

  /**
   * Computes the group with maximum degree centrality of the graph passed in
   * the constructor.
   */
  void run() override;

  /**
   * Returns the group with maximum degree centrality.
   */
  std::vector<node> groupMaxDegree();

  /**
   * Returns the score of the group with maximum degree centrality (i.e. the
   * number of nodes outside the group that can be reached in one hop from at
   * least one node in the group).
   */
  count getScore();

protected:
  Graph G;
  const count k;
  const bool countGroupNodes;
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
  void updateGroup();
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
  groupScore = countGroupNodes ? k : 0;
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
