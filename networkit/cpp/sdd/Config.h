/*
 * Config.h
 *
 *  Created on: May 04, 2014
 *      Author: dhoske
 */

#ifndef SDD_CONFIG_H_
#define SDD_CONFIG_H_

#include <random>
#include <limits>
#include <string>

// Only profile under windows (when callgrind is available)
//#ifndef _WIN32
//#include <valgrind/callgrind.h>
//#else
//#define CALLGRIND_START_INSTRUMENTATION
//#define CALLGRIND_STOP_INSTRUMENTATION
//#define CALLGRIND_DUMP_STATS_AT(x)
//#endif

#include "../algebraic/Vector.h"
#include "../algebraic/Matrix.h"
#include "../graph/Dijkstra.h"
#include "../graph/IncompleteDijkstra.h"

namespace NetworKit {
namespace SDD {


/**
 * @defgroup sdd SDD
 * Solver for linear systems in SDD matrices.
 *
 * This module implements the nearly-linear time solver for SDD-systems
 * described by Kelner et al. in <a href="http://arxiv.org/abs/1301.6628">&ldquo;A
 * Simple, Combinatorial Algorithm for Solving SDD Systems in Nearly-Linear Time&rdquo;</a>.
 *
 * To solve a SDD-system, you should call @link solveSDD<TCycle, TFlow, STAlgo>() @endlink:
 * @code
 *   Matrix A;
 *   Vector b;
 *   SolverStatus status;
 *   status.desired_residual = 1e-3;
 *   Vector x = SDD::solveSDD<UniformCycleDistribution, TrivialFlow, minDistanceST>(A, b, status);
 * @endcode
 * This returns a vector @a x with \f$\Vert Ax - b\Vert_2 / \Vert b \Vert_2 < status.desired\_residual\f$
 * if @a status.converged is true.
 *
 * The @ref solveSDD routine first converts the system into a Laplacian system which
 * is then solved with @ref solveLaplacian. Since this increases the size of the matrix by a
 * factor of 2, you should use @ref solveLaplacian directly if you have
 * a Laplacian system.
 *
 * The solver is parametrised by various algorithms and data structures:
 *   - @a TCycle, the distribution on basis cycles, see @ref sdd-cycle
 *   - @a TFlow, the data structure to store tree flows, see @ref sdd-flow
 *   - @a STAlgo, the spanning tree to use as cycle basis, see @ref sdd-st
 *   - @a status, run-time settings such as the desired residual and
 *     the maximum number of iterations, see @ref sdd-settings
 *
 * @ref sdd-settings also contains several other compile-time settings and
 * @ref sdd-misc contains helper classes and functions.
 *
 * See the thesis by Daniel Hoske for experimental results
 * about this implementation.
 * @todo Add link to thesis
 *
 * @warning We define the stretch
 * of an edge \f$uv\f$ with respect to a spanning tree \f$T\f$ as
 *   \f[st(uv) = \frac{\sum_{e \in P_T(u, v)} 1/w(e)}{1/w(uv)}\f]
 * where \f$P_T(u, v)\f$ denotes the unique path in \f$T\f$ from \f$u\f$
 * to \f$v\f$ and \f$w(e)\f$ is the weight of the edge \f$e\f$.
 */

/** @defgroup sdd-misc Misc
 *  @ingroup sdd
 *
 * Helper classes and functions used in the SDD solvers.
 *
 * This module contains some helper algorithms used in the SDD and Laplacian solvers.
 * In particular:
 *  - @ref LCA is a lowest common ancestor implementation with \f$O(n\log(n))\f$ time precomputation
 * and \f$O(1)\f$ time queries
 *  - @ref RandomTreeGenerator builds random rooted trees using Pr��fer sequences
 *  */

/**
 * @defgroup sdd-settings Settings
 * @ingroup sdd
 * Settings of the SDD-solver.
 *
 * @link SolverStatus @endlink allows you to set configuration options
 * of the solvers and get status information back from them.
 *
 * You can also set various compile-time settings for the solvers.
 */

/** @addtogroup sdd-settings
 *  @{ */

/** Floating point epsilon to use in comparisons. */
constexpr double EPSILON = 1e-9;

/** How many loop iterations to unroll in optimized flow implementations. */
constexpr count UNROLL_LOOP = 4;

/** Use O(1) expected time discrete distribution algorithm? */
#define DISCRETE_O1 0

/** SSSP algorithm to use in the low-stretch spanning tree calculation
 *  when all distances are required. */
using ConeSSSP = Dijkstra;
/** SSSP algorithm to use in the low-stretch spanning tree calculation
 *  when growing balls and not all distances may be required. */
using IncompleteConeSSSP = IncompleteDijkstra;

/** At how many nodes do we switch from complex low-stretch algorithm to a simpler ST-algorithm? */
constexpr count LOW_STRETCH_SMALL = 2;
static_assert(LOW_STRETCH_SMALL == 2, "Switch from Low-stretch to simpler currently only implemented for 2 nodes");

/** Type of a PRNG seed. */
using RandomEngine = std::mt19937;
using SeedType = RandomEngine::result_type;
constexpr SeedType DEFAULT_SEED = RandomEngine::default_seed;

/**
 * Configuration options (input) and status (output) of the SDD solver or the Laplacian solver.
 */
struct SolverStatus {
  /** *Input:* Desired relative residual for the linear system $Lx = y$. */
  double desired_residual = EPSILON;
  /** *Output:* Achieved residual. */
  double actual_residual;

  /** *Input:* Precondition the solver by scaling the
   * spanning tree as in the FullSolver. */
  bool precondition = false;

  /** *Input:* How often is the residual updated.
   *  This is measured in steps of the main loop and
   *  0 stands for an automatic polic. */
  int residual_update_freq = 0;

  /** *Output:* Sum of stretches of the used spanning trees (one for each component).  */
  double stretch;

  /** *Input:* Seed used for the random cycle selection. */
  SeedType seed = DEFAULT_SEED;

  /** *Output:* Has the solver converged?  */
  bool converged;

  /** *Input:* Maximum number of iterations. */
  uint64_t max_iters = std::numeric_limits<uint64_t>::max();
  /**
   * *Output:* Number of computed iterations with the bound in the paper.
   *
   * @remark This is not directly comparable to @link niters @endlink since
   * in the paper we iterate until \f$\Vert x - L^+b \Vert_L < \epsilon\f$,
   * while we run the algorithm until \f$ \Vert Lx - b \Vert_2 / \Vert b \Vert_2 < \epsilon \f$
   */
  uint64_t computed_niters;
  /** *Output:* Actual number of iterations. */
  uint64_t niters;

  /** *Output:* Time for initialising the data structures (in ms). */
  uint64_t time_init;
  /** *Output:* Time for the main loop (in ms). */
  uint64_t time_main;

  /** Debugging fields @{ */

  /** Improvement of residual vs length of the cycle */
  bool debug_enable_stretch_vs_energy = false;
  /** Cycle stretch */
  std::vector<double> debug_stretch;
  /** Energy improvement (for every considered cycle!) */
  std::vector<double> debug_stretch_energy;
  /** Sparsity for every considered cycle */
  std::vector<double> debug_stretch_time;

  /** Trace the development of the residual? */
  bool debug_trace = false;
  /** Trace of residual (only makes sense if the graph is connected) */
  std::vector<double> debug_residuals;
  /** Trace of energy (only makes sense if the graph is connected) */
  std::vector<double> debug_energy;
  /** Trace of iteration count */
  std::vector<count> debug_iters;
  /** Trace of times since start of main loop */
  std::vector<count> debug_times;
  /** Maximum time to execute (in ms) */
  count max_time = count(-1);

  /** Start with given vector?
   *  Requires O(m) additional time at the start.
   *  @bug: Only works if input graph is connected.
   */
  bool debug_enable_initial = false;
  Vector debug_initial;

  /** Trace the computed solutions? */
  bool debug_enable_vectors = false;
  std::vector<Vector> debug_vectors;

  /** @} */
};

/** Prototype of a solver function function. */
using SolverFunc = Vector(const Graph& G, const Vector& b, SolverStatus&);

/****** General purpose helper functions. ******/

/**
 * Computes the relative residual \f$\Vert Lx - b\Vert_2 / \Vert b\Vert_2\f$.
 */
double residual(const Matrix& L, const Vector& x, const Vector& b);

/**
 * Sets the conditions to use during benchmarking.
 * There is currently no way to revert these conditions.
 */
void setBenchmarkConditions();

/** @} */

/** @ingroup sdd-misc
 *  @{ */

/**
 * Computes the sum of elements in the vector @a b.
 */
double vectorSum(const Vector& b);

/**
 * Computes the new names when indexes in @a I are mapped to \f$[0, |I|-1]\f$.
 */
std::unordered_map<node, node> inducedNames(const std::vector<node>& V);

/**
 * Returns the connected components of @a G. This method does not provide
 * the original names of the nodes in the components, i.e. you cannot restore
 * @a G from the components.
 */
std::vector<Graph> getComponents(const Graph& G);

/**
 * Returns the subgraph of @a G induced by the vertices in @a V.
 * (The nodes are taken in the specified order and weights are preserved,
 * unlike @link Subgraph::fromNodes() @endlink).
 */
Graph inducedSubgraph(const Graph& G, const std::vector<node>& V);

/**
 * Returns the subvector of @a v induced by the vertices in @a I.
 */
Vector inducedVector(const Vector& v, const std::vector<node>& I);

/**
 * Returns a random vector whose elements sum to zero
 * in every component of @a G.
 */
Vector randZeroSum(const Graph& G, size_t seed);

/**
 * Returns a graph with inverted weights. This turns a graph
 * containing conductances as weights into a graph containing
 * resistances as weights.
 */
Graph invertGraph(const Graph& G);

/**
 * Reads a graph from @a path and automatically which reader
 * to call based on the extension of @a path.
 */
Graph readGraph(const std::string& path);

/**
 * Generate unweighted 2D-grid of size \f$n \times n\f$.
 * @todo Move to graph generators.
 */
Graph gen2DGrid(count n);

/** @} */

}
}

#endif
