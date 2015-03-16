/*
 * RandomEdgePartitioner.h
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef RANDOMEDGEPARTITIONER_H_
#define RANDOMEDGEPARTITIONER_H_

#include <utility>

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 * Partitions the set of edges of an given graph into two separate edge-sets
 * that get encapsulated into two corresponding graphs. This is done by randomly
 * removing edges from the given graph until a given percentage of edges have been removed.
 */
class RandomEdgePartitioner {
private:
  const Graph& G;

public:
  /**
   * @param G The graph whose edges to partition
   */
  RandomEdgePartitioner(const Graph& G);
  
  /**
   * Randomly removes edges until the given percentage of total edges has been removed.
   * @param percentage Percentage of edges to remove from the graph
   * @return a pair of new graphs where the first graph is the remaining
   * graph and the second graph consists of all removed edges
   */
  std::pair<Graph, Graph> partition(double percentage);

};

} // namespace NetworKit

#endif /* RANDOMEDGEPARTITIONER_H_ */