/*
 * NeighborhoodUtility.h
 *
 *  Created on: 06.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef NEIGHBORHOODUTILITY_H_
#define NEIGHBORHOODUTILITY_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Provides basic operations on neighborhoods in a given graph.
 */
class NeighborhoodUtility {
private:
  /**
   * Sorts and returns the neighborhoods of the nodes @a u and @a v in the given graph @a G.
   * @param G Graph to obtain neighborhoods from
   * @param u First node in the graph
   * @param v Second node in the graph
   * @return the sorted nodes regarding the neighborhoods of the nodes @a u and @a v in the given graph @a G.
   */
  static std::pair<std::vector<node>, std::vector<node>> getSortedNeighborhoods(const Graph& G, node u, node v);

public:
  /**
   * Returns the union of the neighboorhoods of @a u and @a v.
   * @param G Graph to obtain neighbors-union from
   * @param u First node
   * @param v Second node
   * @return a vector containing all the nodes in the neighboorhood-union of @a u and @a v
   */
  static std::vector<node> getNeighborsUnion(const Graph& G, node u, node v);

  /**
   * Returns a vector containing the node-ids of all common neighbors of @a u and @a v.
   * @param G Graph to obtain common neighbors from
   * @param u First node
   * @param v Second node
   * @return a vector containing the node-ids of all common neighbors of @a u and @a v
   */
  static std::vector<node> getCommonNeighbors(const Graph& G, node u, node v);

};

} // namespace NetworKit

#endif /* NEIGHBORHOODUTILITY_H_ */