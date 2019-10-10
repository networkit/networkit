/*
 * UDegreeIndex.h
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef NETWORKIT_LINKPREDICTION_U_DEGREE_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_U_DEGREE_INDEX_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>

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

#endif // NETWORKIT_LINKPREDICTION_U_DEGREE_INDEX_HPP_
