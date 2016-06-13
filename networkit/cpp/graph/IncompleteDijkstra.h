/*
 * IncompleteDijkstra.h
 *
 *  Created on: 15.07.2014
 *      Author: dhoske
 */

#ifndef INCOMPLETEDIJKSTRA_H_
#define INCOMPLETEDIJKSTRA_H_

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <utility>

#include "Graph.h"
#include "IncompleteSSSP.h"
#include "../auxiliary/PrioQueue.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Implementation of @a IncompleteSSSP using a normal
 * Dijkstra with binary heaps.
 */
class IncompleteDijkstra : public IncompleteSSSP {
public:
  /**
   * Creates a IncompleteDijkstra instance from the sources in
   * @a sources (act like a super source) in the graph @a G.
   * The edges in @a G must have nonnegative weight and @a G should
   * not be null.
   *
   * We also consider the nodes in @a explored to not exist
   * if @a explored is not null.
   *
   * @warning We do not copy @a G or @a explored, but store a
   * non-owning pointer to them. Otherwise IncompleteDijkstra would not
   * be more efficient than normal Dijkstra. Thus, @a G and @a explored
   * must exist at least as long as this IncompleteDijkstra instance.
   *
   * @todo This is somewhat ugly, but we do not want introduce a
   * std::shared_ptr<> since @a G and @a explored could well
   * be stack allocated.
   */
  IncompleteDijkstra(const Graph* G, const std::vector<node>& sources,
                     const std::unordered_set<node>* explored = nullptr);

  virtual bool hasNext() override;
  virtual std::pair<node, edgeweight> next() override;

private:
  // discard duplicate elements in pq
  void discardDuplicates();

  // Stored reference to outside data structures
  const Graph* G;
  const std::unordered_set<node>* explored;

  // distances aren't stored in a vector because initialising it may be too expensive
  std::unordered_map<node, edgeweight> dists;
  // TODO: Fix Aux::PrioQueue to work with arbitrary values.
  // and use it instead.
  using PrioValue = std::pair<edgeweight, node>;
  using Prio = std::priority_queue<PrioValue, std::vector<PrioValue>, std::greater<PrioValue>>;
  Prio pq;
};

}

#endif
