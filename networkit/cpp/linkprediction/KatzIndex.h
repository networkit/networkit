/*
 * KatzIndex.h
 *
 *  Created on: 30.01.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef KATZINDEX_H_
#define KATZINDEX_H_

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 * Katz index used for similarity measurements on graphs.
 */
class KatzIndex : public LinkPredictor {
private:
  unsigned int maxPathLength;

  double dampingValue;

  node lastStartNode;

  std::unordered_map<node, double> lastScores;

public:
  /**
   * @param maxPathLength limit for the length of the paths to consider
   */
  KatzIndex(const Graph& G, unsigned int maxPathLength, double dampingValue);
  
  double run(node u, node v);

};

} // namespace NetworKit

#endif /* KATZINDEX_H_ */