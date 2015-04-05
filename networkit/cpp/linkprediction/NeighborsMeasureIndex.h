/*
 * NeighborsMeasureIndex.h
 *
 *  Created on: 05.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef NEIGHBORSMEASUREINDEX_H_
#define NEIGHBORSMEASUREINDEX_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * ...
 */
class NeighborsMeasureIndex : public LinkPredictor {
private:

  /**
   * Returns the number of connections between neighbors of @a u and @a v.
   * @param u First node
   * @param v Second node
   * @return the number of connections between neighbors of @a u and @a v
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;
  
};

} // namespace NetworKit

#endif /* NEIGHBORSMEASUREINDEX_H_ */