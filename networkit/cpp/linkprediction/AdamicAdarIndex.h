/*
 * AdamicAdarIndex.h
 *
 *  Created on: 25.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef ADAMICADAR_H_
#define ADAMICADAR_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Implementation of the Adamic/Adar Index.
 * The index sums up the reciprocals of the logarithm of the degree of all
 * common neighbors of u and v.
 */
class AdamicAdarIndex : public LinkPredictor {
private:
  /**
   * Returns the Adamic/Adar Index of the given node-pair (@a u, @a v).
   * @param u First node
   * @param v Second node
   * @return the Adamic/Adar Index of the given node-pair (@a u, @a v)
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif /* ADAMICADAR_H_ */