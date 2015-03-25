/*
 * AdamicAdarIndex.h
 *
 *  Created on: 25.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef ADAMICADAR_H_
#define ADAMICADAR_H_

#include "LinkPredictor.h"
#include "CommonNeighborsIndex.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 *
 */
class AdamicAdarIndex : public LinkPredictor {
private:
  CommonNeighborsIndex commonNeighborsIndex;

  /**
   *
   * @param u First node
   * @param v Second node
   * @return 
   */
  double runImpl(node u, node v) override;

public:
  AdamicAdarIndex();

  /**
   *
   * @param G The graph to work on
   */
  explicit AdamicAdarIndex(const Graph& G);
};

} // namespace NetworKit

#endif /* ADAMICADAR_H_ */