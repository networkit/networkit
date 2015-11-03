/*
 * VDegreeIndex.h
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef VDEGREEINDEX_H_
#define VDEGREEINDEX_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Index that simply returns the degree of the second given node.
 */
class VDegreeIndex : public LinkPredictor {
private:
  /**
   * Returns the degree of the second node provided, namely @a v.
   * @param u First node
   * @param v Second node
   * @return the degree of the second node provided, namely @a v
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;
  
};

} // namespace NetworKit

#endif /* VDEGREEINDEX_H_ */