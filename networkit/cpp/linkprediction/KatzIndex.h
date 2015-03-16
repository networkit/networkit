/*
 * KatzIndex.h
 *
 *  Created on: 30.01.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef KATZINDEX_H_
#define KATZINDEX_H_

#include <unordered_map>

#include "LinkPredictor.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 * Katz index assigns a pair of nodes a similarity score
 * that is based on the sum of the weighted number of paths of length l
 * where l is smaller than a given limit.
 */
class KatzIndex : public LinkPredictor {
private:
  unsigned int maxPathLength; //!< Maximal length of the paths to consider

  double dampingValue; //!< Damping factor in (0,1) used to exponentially damp every addend of the sum

  node lastStartNode; //!< Node at which the algorithm started at the last call

  // TODO: Maybe replace with vector/array of size G.upperNodeIdBound? Constant access but uses more memory.
  std::unordered_map<node, double> lastScores;

  // Helper method used to access the score for a given node-pair. Checks which of the given nodes
  // was used as the starting node and uses the other node to access the last scores generated.
  // Defaults to 0.0 if no score could be found.
  double getScore(node u, node v) const;

  bool executed = false;

public:
  /**
   * @param G The graph to operate on
   * @param maxPathLength Maximal length of the paths to consider
   * @param dampingValue Used to exponentially damp every addend of the sum. Should be in (0, 1)
   */
  KatzIndex(const Graph& G, count maxPathLength = 3, double dampingValue = 0.9);
  
  /**
   * Returns the similarity score for the given node-pair based on the
   * Katz index specified during construction. The algorithm considers all
   * paths starting at the node with the smaller degree except the algorithm
   * started at the other node at the last call.
   * @param u First node
   * @param v Second node
   * @return the similarity score of the given node-pair calculated by the specified Katz index
   */
  double run(node u, node v);

};

} // namespace NetworKit

#endif /* KATZINDEX_H_ */