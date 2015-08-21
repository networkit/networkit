/*
 * CommonNeighborsIndex.h
 *
 *  Created on: 06.12.2014
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef COMMONNEIGHBORSINDEX_H_
#define COMMONNEIGHBORSINDEX_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * The CommonNeighborsIndex calculates the number of common
 * neighbors of a node-pair in a given graph.
 */
class CommonNeighborsIndex : public LinkPredictor {
private:
  /**
   * Returns the number of common neighbors of the given nodes @a u and @a v.
   * @param u First node
   * @param v Second node
   * @return the number of common neighbors of @a u and @a v
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif /* COMMONNEIGHBORSINDEX_H_ */