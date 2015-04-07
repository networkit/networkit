/*
 * TotalNeighborsIndex.h
 *
 *  Created on: 05.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef TOTALNEIGHBORSINDEX_H_
#define TOTALNEIGHBORSINDEX_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Implementation of the Total Neighbors Index.
 * This index is also known as Total Friends Index and returns
 * the number of nodes in the neighborhood-union of u and v.
 */
class TotalNeighborsIndex : public LinkPredictor {
private:
  /**
   * Returns the number of total union-neighbors for the given node-pair (@a u, @a v).
   * @param u First node
   * @param v Second node
   * @return the number of total union-neighbors for the given node-pair (@a u, @a v)
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif /* TOTALNEIGHBORSINDEX_H_ */