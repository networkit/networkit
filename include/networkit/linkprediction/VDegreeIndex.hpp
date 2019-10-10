/*
 * VDegreeIndex.h
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef NETWORKIT_LINKPREDICTION_V_DEGREE_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_V_DEGREE_INDEX_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>

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

#endif // NETWORKIT_LINKPREDICTION_V_DEGREE_INDEX_HPP_
