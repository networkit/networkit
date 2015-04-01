/*
 * TrainingGraphGenerator.h
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef TRAININGGRAPHGENERATOR_H_
#define TRAININGGRAPHGENERATOR_H_

#include <utility>

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Partitions the set of edges of an given graph into two separate edge-sets
 * that get encapsulated into two corresponding graphs. This is done by randomly
 * removing edges from the given graph until a given percentage of edges have been removed.
 */
class TrainingGraphGenerator {
public:
  /**
   * Randomly removes edges until the given percentage of total edges has been removed.
   * Removed edges will be added to a new graph.
   *
   * @param percentage Percentage of edges to remove from the graph
   * @return a pair of new graphs where the first graph is the remaining
   * graph and the second graph consists of all removed edges
   */
  static Graph byPercentage(const Graph& G, double trainPercentage);

  /**
   * Randomly removes edges until the given count of total edges has been removed.
   * Removed edges will be added to a new graph.
   *
   * @param numEdges Number of edges to remove from the graph
   * @return a pair of new graphs where the first graph is the remaining
   * graph and the second graph consists of all removed edges
   */
  static Graph byCount(const Graph& G, count numTrainEdges);

};

} // namespace NetworKit

#endif /* TRAININGGRAPHGENERATOR_H_ */