/*
 * CycleDistribution.h
 *
 *  Created on: Apr 26, 2014
 *      Author: dhoske
 */

#ifndef CYCLE_DISTRIBUTION_H_
#define CYCLE_DISTRIBUTION_H_

#include <random>

#include "../graph/Graph.h"
#include "Config.h"
#include "RootedTree.h"

namespace NetworKit {
namespace SDD {


/**
 * @defgroup sdd-cycle Cycle distribution
 * @ingroup sdd
 * Distributions on the basis cycles induced by a spanning tree of a graph.
 *
 * The classes derived from @ref CycleDistribution return
 * random basis cycle induced by a spanning tree of a graph.
 * @ref UniformCycleDistribution choses cycles uniformly at random,
 * while @ref StretchCycleDistribution weighs them by their stretch.
 *
 *  @{ */

/**
 * Provides a distribution on the basis cycles given by
 * a spanning tree of a graph.
 *
 * Note: randomCycle() occurs in the inner loop of the algorithm and
 * should be inlined.
 */
class CycleDistribution {
public:
  /**
   * Creates a cycle distribution on the cycle basis given
   * by the spanning tree @a T on the graph @a G.
   *
   * We must live with increased coupling for reasons of efficiency:
   * The vector @a off_tree_edges should contain all of the off-tree-edges
   * in @a G with respect to to @a G. We could compute this on our own,
   * but we want to have a consistent ordering of the off-tree-edges
   * in order to be able to index off-tree-edges by their index in @a off_tree_edges
   * and avoid hashing. (This optimization improves the running time by ~40%.)
   * We also provide the stretches of the edges from the outside, so they aren't recomputed.
   */
  CycleDistribution(const Graph& G, const RootedTree& T,
                    const std::vector<Edge>& off_tree_edges, const std::vector<edgeweight>& stretches,
                    SeedType seed);

  /**
   * Returns a random basis cycles (off-tree-edge) in the graph represented by
   * its index in the vector of off-tree-edges given in @ref CycleDistribution.
   */
  virtual index randomCycle() = 0;

  /**
   * Returns the (unnormalized) probability of chosing edge @a e.
   */
  virtual double getProb(const Edge& e) const;

  virtual ~CycleDistribution() = default;

  /* Neither copyable nor moveable */
  CycleDistribution(const CycleDistribution&) = delete;
  CycleDistribution& operator=(const CycleDistribution&) = delete;
  CycleDistribution(CycleDistribution&&) = delete;
  CycleDistribution& operator=(CycleDistribution&&) = delete;

protected:
  RandomEngine engine;
  std::vector<Edge> off_tree_edges;
};

/**
 * Cycle distribution that choses each cycle uniformly at random.
 */
class UniformCycleDistribution : public CycleDistribution {
public:
  /** See CycleDistribution::CycleDistribution() */
  UniformCycleDistribution(const Graph& G, const RootedTree& T, const std::vector<Edge>& off_tree_edges, const std::vector<edgeweight>& stretches, SeedType seed);

  virtual index randomCycle() override {
    return random_edge(engine);
  }

private:
  std::uniform_int_distribution<index> random_edge;
};

/**
 * Cycle distribution that weighs each cycle by its stretch.
 */
class StretchCycleDistribution : public CycleDistribution {
public:
  /** See CycleDistribution::CycleDistribution(). */
  StretchCycleDistribution(const Graph& G, const RootedTree& T, const std::vector<Edge>& off_tree_edges, const std::vector<edgeweight>& stretches, SeedType seed);

  virtual index randomCycle() override {
#if(DISCRETE_O1)
    auto random_index = std::uniform_int_distribution<size_t>(0, probs.size() - 1);
    auto random_prob  = std::uniform_real_distribution<double>(0.0, max_prob);
    index out = 0;
    while (true) {
      out = random_index(engine);
      if (random_prob(engine) < probs[out]) {
        break;
      }
    }
    return out;
#else
    return random_edge(engine);
#endif
  }

  virtual double getProb(const Edge& e) const override;

private:
  std::vector<double> probs;
  std::discrete_distribution<index> random_edge;
#if(DISCRETE_O1)
  double max_prob;
#endif
};

/** @} */

}
}

#endif
