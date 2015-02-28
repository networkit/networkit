/*
 * RandomEdgeRemover.h
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef RANDOMEDGEREMOVER_H_
#define RANDOMEDGEREMOVER_H_

#include <utility>

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 * Returns a new graph where a given percentage of edges has been randomly removed.
 */
 // TODO: Maybe rename to RandomGraphPartitioner
class RandomEdgeRemover {
private:
  const Graph& G;

public:
  /**
   * @param G Graph to remove edges from
   */
  RandomEdgeRemover(const Graph& G);
  
  /**
   * Randomly removes edges until the given percentage of
   * total edges has been removed.
   * @return a pair of new graphs where the first graph is the remaining
   * graph and the second graph consists of all removed edges
   */
  std::pair<Graph, Graph> remove(double percentage);

};

} // namespace NetworKit

#endif /* RANDOMEDGEREMOVER_H_ */