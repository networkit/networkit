/*
 * TotalNeighborsIndex.h
 *
 *  Created on: 05.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef TOTALNEIGHBORSINDEX_H_
#define TOTALNEIGHBORSINDEX_H_

#include "LinkPredictor.h"
#include "JaccardIndex.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * ...
 */
class TotalNeighborsIndex : public LinkPredictor {
private:
  JaccardIndex jaccardIndex; //!< Used to get the neighborhood-union of u and v

  /**
   * Returns the number of total union-neighbors for the given node-pair (@a u, @a v).
   * @param u First node
   * @param v Second node
   * @return the number of total union-neighbors for the given node-pair (@a u, @a v)
   */
  double runImpl(node u, node v) override;

public:
  TotalNeighborsIndex();

  explicit TotalNeighborsIndex(const Graph& G);
  
  void setGraph(const Graph& newGraph) override;
};

} // namespace NetworKit

#endif /* TOTALNEIGHBORSINDEX_H_ */