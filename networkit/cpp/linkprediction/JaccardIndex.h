/*
 * JaccardIndex.h
 *
 *  Created on: 23.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef JACCARDINDEX_H_
#define JACCARDINDEX_H_

#include "LinkPredictor.h"
#include "CommonNeighborsIndex.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Implementation of the Jaccard index which normalizes the Common Neighbors Index
 * by dividing through the number of nodes in the neighboorhood-union.
 */
class JaccardIndex : public LinkPredictor {
private:
  CommonNeighborsIndex commonNeighborsIndex; //!< Used to get the number of common neighbors of u and v

  /**
   * Returns the Jaccard index for the given node-pair (@a u, @a v).
   * @param u First node
   * @param v Second node
   * @return the Jaccard index for the given node-pair (@a u, @a v)
   */
  double runImpl(node u, node v) override;

public:
  JaccardIndex();

  /**
   *
   * @param G The graph to work on
   */
  explicit JaccardIndex(const Graph& G);

  /**
   * Returns the union of the neighboorhoods of @a u and @a v.
   * @param u First node
   * @param v Second node
   * @return a vector containing all the nodes in the neighboorhood-union of u and v
   */
  std::vector<node> getNeighborsUnion(node u, node v) const;
};

} // namespace NetworKit

#endif /* JACCARDINDEX_H_ */