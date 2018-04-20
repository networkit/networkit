/*
 * GroupDegree.h
 *
 *  Created on: 20.04.2018
 *      Author: Eugenio Angriman
 */

#ifndef GROUPDEGREE_H_
#define GROUPDEGREE_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

class GroupDegree : public Algorithm {
public:
  GroupDegree(const Graph &G, count k = 1);
  void run() override;
  std::vector<node> groupMaxDegree();

protected:
  Graph G;
  const count k;
};
} // namespace NetworKit

#endif
