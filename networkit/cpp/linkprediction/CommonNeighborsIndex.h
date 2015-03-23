/*
 * CommonNeighborsIndex.h
 *
 *  Created on: 06.12.2014
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef COMMONNEIGHBORSINDEX_H_
#define COMMONNEIGHBORSINDEX_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * The CommonNeighborsIndex calculates the number of common
 * neighbors of a node-pair in a given graph.
 */
class CommonNeighborsIndex : public LinkPredictor {
private:
  /**
   * Returns the number of common neighbors of the given nodes @a u and @a v.
   * @param u First node
   * @param v Second node
   * @return the number of common neighbors of u and v
   */
  double runImpl(node u, node v) override;

public:
  CommonNeighborsIndex();

  /**
   *
   * @param G The graph to work on
   */
  explicit CommonNeighborsIndex(const Graph& G);

  /**
   * Returns the common neighbors of the given nodes @a u and @a v.
   * 
   *
   */
  std::vector<node> getCommonNeighbors(node u, node v) const;

};

} // namespace NetworKit

#endif /* COMMONNEIGHBORSINDEX_H_ */