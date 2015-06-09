/*
 * LCA.h
 *
 *  Created on: 19.07.2014
 *      Author: dhoske
 */

#ifndef LCA_H_
#define LCA_H_

#include <random>

#include "RootedTree.h"

namespace NetworKit {
namespace SDD {

/** @ingroup sdd-misc
 *  @{ */

/**
 * Data structure for answering LCA queries in a tree with
 * \f$O(n\log(n))\f$ time precomputation and \f$O(1)\f$ time
 * queries using sparse tables.
 */
class LCA {
public:
  /**
   * Builds the LCA data structure on @a T in \f$O(n\log(n))\f$ time.
   */
  explicit LCA(const RootedTree& T);

  /**
   * Answer an LCA query between @a u and @a v in \f$O(n)\f$ time.
   */
  node query(node u, node v) const;

  /* Copyable and moveable */
  LCA(const LCA&) = default;
  LCA(LCA&&) = default;
  LCA& operator=(const LCA&) = default;
  LCA& operator=(LCA&&) = default;

private:
  // Number of nodes
  count n;
  // First occurence of a node in Eulerian tour
  std::vector<index> first_occ;
  // Depths of nodes
  std::vector<count> depths;
  // RMQ data structure (for each length that is a power of two we store the precomputed RMQs)
  std::vector<std::vector<node>> rmq;
};

/** @} */

}
}

#endif
