/*
 * Flow.h
 *
 *  Created on: May 03, 2014
 *      Author: dhoske
 */

#ifndef FLOW_H_
#define FLOW_H_

#include <climits>
#include "../algebraic/Vector.h"
#include "Config.h"
#include "RootedTree.h"
#include "LCA.h"

//#ifndef NOVECTORIZE
//#include "vectorclass.h"
//#endif

namespace NetworKit {
namespace SDD {

/** @defgroup sdd-flow Flow
 *  @ingroup sdd
 *
 *  Data structures for storing and updating flows in a tree.
 *
 *  The classes derived from @ref Flow store flows in a tree.
 *
 *  You can add flow and get the potential of each node in the tree (weights
 *  are interpreted as conductances). @ref TrivialFlow
 *  implements all these operations in \f$O(n)\f$ time in the worst case,
 *  while @link LogFlow @endlink only needs \f$O(\log(n))\f$ time for each operation.
 *  @ref LCAFlow is an improvement of @ref TrivialFlow that only walks
 *  to the LCA of the two nodes, i.e. it will be faster than @ref TrivialFlow
 *  but still needs linear time for each operation in the worst case.
 *  @{ */

/**
 * Computes the unique flow with outflows given by @a b that is non-zero only on
 * the edges of the spanning tree @a T. The entries of @a b must sum to zero.
 */
RootedTree initialFlow(const RootedTree& T, const Vector& b);

/**
 * Data structure for storing a flow (currents) in a rooted tree
 * where the edge weights are interpreted as conductances.
 *
 * Note: \ref addFlow(), \ref getPotential() and \ref primalSolution() are
 * in the inner loop of the  Laplace solver, so they should be inlined.
 * (The performance difference if they not inlined is ~5%.)
 */
class Flow {
public:
  /**
   * Constructs a new flow with the topology given by @a T. The edge
   * weights of @a T are interpreted as conductances.
   */
  explicit Flow(const RootedTree& T);

  /**
   * @copydetails Flow(const RootedTree&)
   * Also sets the initial flow to the flow given by the outflows in
   * @a b at each node. The entries in @a b must sum to zero.
   */
  Flow(const RootedTree& T, const Vector& b);

  /**
   * Returns the primal solution (the voltages) represented by the
   * electrical flow with respect to the root of the tree. This requires
   * \f$|T|\f$ getPotential()-operations if you use this generic implementation.
   */
  virtual Vector primalSolution() const;

  /**
   * Returns the voltage drop along the path from @a u to @a v.
   */
  virtual edgeweight getPotential(node u, node v) const = 0;

  /**
   * Adds @a amount flow to the path from @a u to @a v.
   */
  virtual void addFlow(node u, node v, edgeweight amount) = 0;

  /**
   * Returns the time required for updating and querying node @a u
   * in some implementation-dependent measure.
   */
  virtual count getTime(node u) const {
    return 0;
  };

  /**
   * Scales the resistances in this flow by @a factor.
   */
  virtual void scaleResistances(edgeweight factor);

  virtual ~Flow() = default;

  /* Copyable and moveable */
  Flow(const Flow&) = default;
  Flow& operator=(const Flow&) = default;
  Flow(Flow&&) = default;
  Flow& operator=(Flow&&) = default;

protected:
  /* Original spanning tree with reciprocal of edge weights (= resistances) on edges */
  RootedTree Tresistances;
};

/**
 * Trivial flow implementation that needs \f$O(n)\f$ worst-case time
 * for each operation.
 */
class TrivialFlow : public Flow {
public:
  /**
   * See Flow::Flow
   */
  TrivialFlow(const RootedTree& T, const Vector& b);

  virtual void addFlow(node u, node v, edgeweight amount) override {
    addFlow(u,  amount);
    addFlow(v, -amount);
  }

  virtual edgeweight getPotential(node u, node v) const override {
    return getPotential(u) - getPotential(v);
  }

  virtual Vector primalSolution() const override;

protected:
  // Get potential with respect to the root
  inline edgeweight getPotential(node u) const;

  // Add flow along u->root path
  inline void addFlow(node u, edgeweight amount);

  /* Tree that stores the flow values at its edges. */
  RootedTree flowTree;
};

/**
 * Flow implementation using LCA that needs \f$O(n)\f$ worst-case time
 * for each operation, but is usually faster than @ref TrivialFlow since
 * each operation only walks to the LCA of the two nodes.
 *
 * May be used in parallel if the two cycles that are modified are disjoint.
 */
class LCAFlow : public TrivialFlow {
public:
  /**
   * See Flow::Flow
   */
  LCAFlow(const RootedTree& T, const Vector& b);

  virtual void addFlow(node u, node v, edgeweight amount) override {
    node lca_uv = lca.query(u, v);
    addFlowUntil(u, lca_uv, amount);
    addFlowUntil(v, lca_uv, -amount);
  }

  virtual edgeweight getPotential(node u, node v) const override {
    node lca_uv = lca.query(u, v);
    return getPotentialUntil(u, lca_uv) - getPotentialUntil(v, lca_uv);
  }

private:
  // Get potential with respect to until
  inline edgeweight getPotentialUntil(node u, node until) const;

  // Add flow along u->until path
  inline void addFlowUntil(node u, node until, edgeweight amount);

  // LCA data structure on T
  LCA lca;
};

/**
 * Improved flow implementation that needs \f$O(\log(n))\f$
 * worst-case time for each operation.
 */
class LogFlow : public Flow {
public:
  /**
   * See Flow::Flow
   */
  LogFlow(const RootedTree& T);

  /**
   * See Flow::Flow
   */
  LogFlow(const RootedTree& T, const Vector& b);

  index getTime(index) const override;

  virtual void scaleResistances(edgeweight factor) override;

  virtual edgeweight getPotential(node u, node v) const override {
    return getPotential(u) - getPotential(v);
  }

  virtual void addFlow(node u, node v, edgeweight amount) override {
    addFlow(u, amount);
    addFlow(v, -amount);
  }

private:
  // Get potential with respect to the root
  inline edgeweight getPotential(node u) const;

  // Add flow along u->root path
  inline void addFlow(node u, edgeweight amount);

  /* State of the data structure (pair of extension and drop values) */
  std::vector<edgeweight> state;

  /* Update vector for each node -> addFlow of u by alpha does state += alpha * update(u) */
  std::vector<std::vector<index>>      update_index;
  std::vector<std::vector<edgeweight>> update_weight;
  /* Query vector for each node -> getPotential of u returns dot(state, query[u]) */
  std::vector<std::vector<index>>      query_index;
  std::vector<std::vector<edgeweight>> query_weight;

  /* Utility functions for the tree decomposition. */
  /* Get sizes of all subtrees of the tree where some nodes can additionally be considered leaves. */
  void getSizes(node root, const std::vector<bool>& is_leaf, std::vector<int>& tree_size /* out */) const;
  /* Find splitting point of the tree of size nnodes starting with the edge e into two almost equal subtrees. */
  node findSplit(const std::vector<int>& tree_size, Edge e, int nnodes) const;
  /* Get heights of nodes (measured in resistance on root-split path) */
  std::vector<edgeweight> getHeights(const std::vector<bool>& is_leaf, node root, node split) const;
};

/***** Inline implementations of operations that are used in the inner loop of the solver *****/
edgeweight TrivialFlow::getPotential(node u) const {
  edgeweight pot = 0;
  flowTree.forEdgesUp([&] (node u, node v, edgeweight w) {
    /* Potential = Current / Conductance */
    pot += w * Tresistances.getWeight(v, u);
  }, u);
  return pot;
}

void TrivialFlow::addFlow(node u, edgeweight amount) {
  flowTree.forEdgesUp([&] (node u, node v, edgeweight weight) {
    flowTree.setWeight(v, u, weight + amount);
  }, u);
}

edgeweight LCAFlow::getPotentialUntil(node u, node until) const {
  edgeweight pot = 0;
  flowTree.forEdgesUpUntil([&] (node u, node v, edgeweight w) {
    /* Potential = Current / Conductance */
    pot += w * Tresistances.getWeight(v, u);
  }, u, until);
  return pot;
}

void LCAFlow::addFlowUntil(node u, node until, edgeweight amount) {
  flowTree.forEdgesUpUntil([&] (node u, node v, edgeweight weight) {
    flowTree.setWeight(v, u, weight + amount);
  }, u, until);
}

edgeweight LogFlow::getPotential(node u) const {
  const auto& qi = query_index[u];
  const auto& qw = query_weight[u];
  assert(qi.size() == qw.size());

//#ifndef NOVECTORIZE
//  // Manually vectorized implementation with Agner's vector classes
//  // (~10% improvement over automatically vectorized implementation)
//  // Still very cache-inefficient.
//  static_assert(UNROLL_LOOP == 4, "Vectorization only supported for 4 elements at a time");
//  Vec4d out(0.0), l, r;
//  Vec4q indexes;
//  for (index i = 0; i < qi.size(); i += 4) {
//    indexes.load(&qi[i]);
//    l = lookup<INT_MAX>(indexes, &state[0]);
//    r.load(&qw[i]);
//    out += l * r;
//  }
//  return horizontal_add(out);
//#else
  edgeweight out = 0;
  //#pragma omp parallel for reduction(+:out)
  for (index i = 0; i < qi.size(); i += UNROLL_LOOP) {
    // This loop should be unrolled when compiled with -O3
    for (index j = 0; j < UNROLL_LOOP; ++j) {
      out += state[qi[i + j]] * qw[i + j];
    }
  }
  return out;
//#endif
}

void LogFlow::addFlow(node u, edgeweight amount) {
  const auto& upi = update_index[u];
  const auto& upw = update_weight[u];

//#pragma omp parallel for
  // Note: vectorization is counterproductive since the
  // access pattern is very irregular, i.e. we only unroll
  // the loop. Still very cache-inefficient.
  for (index i = 0; i < upi.size(); i += UNROLL_LOOP) {
    // This loop should be unrolled when compiled with -O3
    for (index j = 0; j < UNROLL_LOOP; ++j) {
      state[upi[i + j]] += amount * upw[i + j];
    }
  }
}

/** @} */

}
}

#endif
