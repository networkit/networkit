/*
 * BFS.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef REVERSEBFS_H_
#define REVERSEBFS_H_

#include "../distance/SSSP.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup graph
 * The BFS class is used to do a reverse breadth-first search (following
 * the incoming edges of a node) on a Graph from a given source node.
 */
class ReverseBFS : public SSSP {

public:
  /**
   * Constructs the BFS class for @a G and source node @a source.
   *
   * @param G The graph.
   * @param source The source node of the breadth-first search.
   * @param storePaths	store paths and number of paths?
   * @param storeStack	maintain a stack of nodes in decreasing order of
   * distance
   */
  ReverseBFS(const Graph &G, node source, bool storePaths = true,
             bool storeNodesSortedByDistance = false, node target = none);

  /**
   * Reverse Breadth-first search from @a source.
   * @return Vector of unweighted distances from node @a source, i.e. the
   * length (number of edges) of the shortest path from @a source to any other
   * node.
   */
  virtual void run();
};

} /* namespace NetworKit */
#endif /* REVERSEBFS_H_ */
