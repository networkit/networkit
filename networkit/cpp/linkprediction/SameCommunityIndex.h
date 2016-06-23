/*
 * SameCommunityIndex.h
 *
 *  Created on: 07.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef SAMECOMMUNITYINDEX_H_
#define SAMECOMMUNITYINDEX_H_

#include "LinkPredictor.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Index to determine whether two nodes are in the same community.
 */
class SameCommunityIndex : public LinkPredictor {
private:
  Partition communities; //!< The communities of the current graph

  /**
   * Returns 1 if the given nodes @a u and @a v are in the same community, 0 otherwise.
   * @param u First node
   * @param v Second node
   * @return 1 if the given nodes @a u and @a v are in the same community, 0 otherwise
   */
  double runImpl(node u, node v) override;

public:
  SameCommunityIndex();

  explicit SameCommunityIndex(const Graph& G);
  
  virtual void setGraph(const Graph& newGraph) override;

};

} // namespace NetworKit

#endif /* SAMECOMMUNITYINDEX_H_ */