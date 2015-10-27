/*
 * Dijkstra.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include "Graph.h"
#include "SSSP.h"
#include "../auxiliary/PrioQueue.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Dijkstra's SSSP algorithm.
 */
class Dijkstra : public SSSP {

public:

  /**
   * Creates the Dijkstra class for @a G and the source node @a source.
   *
   * @param G The graph.
   * @param source The source node.
   * @param storePaths	store paths and number of paths?
   * @param storeStack
   */
  Dijkstra(const Graph& G, node source, bool storePaths=true, bool storeStack=false);

  /**
   * Creates the Dijkstra class for @a G and the source nodes in @a sources.
   *
   * @param G The graph.
   * @param sources The source nodes.
   */
  Dijkstra(const Graph& G, const std::vector<node>& sources);

  /**
   * Performs the Dijkstra SSSP algorithm on the graph given in the constructor.
   */
  virtual void run(node t = none);

  /**
   * Performs the Dijkstra SSSP algorithm search from @a source until @a target.
   *
   * @param t The target node. The search will stop as soon as its found.
   */
  void runUntil(node t);

};

} /* namespace NetworKit */
#endif /* DIJKSTRA_H_ */
