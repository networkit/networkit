/*
 * ReversedNeighborhoodDistanceIndex.h
 *
 *  Created on: 24.06.2013
 *      Authors: cls, Kolja Esders
 */

#ifndef REVERSEDNEIGHBORHOODDISTANCEINDEX_H_
#define REVERSEDNEIGHBORHOODDISTANCEINDEX_H_

#include "LinkPredictor.h"
#include <math.h>
#include <algorithm>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Assigns a distance value to pairs of nodes according to the
 * overlap of their neighborhoods.
 */
class ReversedNeighborhoodDistanceIndex : public LinkPredictor {
private:
  /**
   * Returns the Neighborhood Distance index for the given node-pair (@a u, @a v).
   * @param u First node
   * @param v Second node
   * @return the Neighborhood Distance index for the given node-pair (@a u, @a v)
   */
  virtual double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;
  
};

} /* namespace NetworKit */
#endif /* REVERSEDNEIGHBORHOODDISTANCEINDEX_H_ */
