/*
 * AdjustedRandIndex.h
 *
 *  Created on: 11.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef NETWORKIT_LINKPREDICTION_ADJUSTED_RAND_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_ADJUSTED_RAND_INDEX_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * AdjustedRandIndex proposed by Hoffman et al. with natural threshold of 0.
 * See http://www.sciencedirect.com/science/article/pii/S0378873315000210 for details.
 */
class AdjustedRandIndex : public LinkPredictor {
private:
  /**
   * Returns the Adjusted Rand Index of the given node-pair (@a u, @a v).
   * @param u First node
   * @param v Second node
   * @return the Adjusted Rand Index of the given node-pair (@a u, @a v)
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_ADJUSTED_RAND_INDEX_HPP_
