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
 *
 */
class JaccardIndex : public LinkPredictor {
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
  JaccardIndex();

  /**
   *
   * @param G The graph to work on
   */
  explicit JaccardIndex(const Graph& G);

  std::vector<node> getNeighborsUnion(node u, node v) const;
};

} // namespace NetworKit

#endif /* JACCARDINDEX_H_ */