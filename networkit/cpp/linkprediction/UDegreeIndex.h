/*
 * UDegreeIndex.h
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef UDEGREEINDEX_H_
#define UDEGREEINDEX_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Index that simply returns the degree of the first given node.
 */
class UDegreeIndex : public LinkPredictor {
private:
  /**
   * Returns the degree of the first node provided, namely @a u.
   * @param u First node
   * @param v Second node
   * @return the degree of the first node provided, namely @a u
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif /* UDEGREEINDEX_H_ */