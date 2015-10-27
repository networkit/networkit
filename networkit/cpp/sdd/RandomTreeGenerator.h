/*
 * RandomTreeGenerator.h
 *
 *  Created on: 19.04.2014
 *      Author: dhoske
 */

#ifndef RANDOMTREEGENERATOR_H_
#define RANDOMTREEGENERATOR_H_

#include <random>
#include "Config.h"

#include "RootedTree.h"

namespace NetworKit {
namespace SDD {

/** @ingroup sdd-misc
 *  @{ */

/**
 * Generator for random ordered rooted trees with edge weights
 * using Prüfer sequences.
 */
class RandomTreeGenerator {
public:
  /**
   * @brief Creates a new tree generator.
   * @param[in] seed seed for the pseudorandom generator
   */
  RandomTreeGenerator(SeedType seed = DEFAULT_SEED) : rand(seed) {
  }

  /**
   * Generates a random rooted tree on @a n vertices with random edge weights
   * in \f$[lower, upper]\f$.
   */
  SDD::RootedTree makeTree(count n, edgeweight lower = 0.0, edgeweight upper = 1.0);

  /* Neither copyable nor moveable */
  RandomTreeGenerator(const RandomTreeGenerator&) = delete;
  RandomTreeGenerator(RandomTreeGenerator&&) = delete;
  RandomTreeGenerator& operator=(const RandomTreeGenerator&) = delete;
  RandomTreeGenerator& operator=(RandomTreeGenerator&&) = delete;

private:
  RandomEngine rand;

  /**
   * Converts the Prüfer sequence @a pruefer to a tree rooted on @a root with
   * @a n nodes and edge weights taken from @a weights.
   */
  RootedTree prueferToTree(count n, node root,
                           const std::vector<node>& pruefer,
                           const std::vector<edgeweight>& weights);
};

/** @} */

}
}

#endif
