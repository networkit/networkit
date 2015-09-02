/*
 * ResourceAllocationIndex.h
 *
 *  Created on: 11.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef RESOURCEALLOCATIONINDEX_H_
#define RESOURCEALLOCATIONINDEX_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Implementation of the ResourceAllocationIndex.
 * The index is similar to Adamic/Adar and sums up the reciprocals of
 * the degree of all common neighbors of u and v.
 */
class ResourceAllocationIndex : public LinkPredictor {
private:
  /**
   * Returns the Resource Allocation Index of the given node-pair (@a u, @a v).
   * @param u First node
   * @param v Second node
   * @return the Resource Allocation Index of the given node-pair (@a u, @a v)
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif /* RESOURCEALLOCATIONINDEX_H_ */