/*
 * LaplaceSolver.h
 *
 *  Created on: May 03, 2014
 *      Author: dhoske
 */

#ifndef LAPLACIAN_SOLVER_H
#define LAPLACIAN_SOLVER_H

#include <vector>
#include <unordered_map>
#include <iomanip>
#include <random>

#include "../graph/Graph.h"
#include "../algebraic/Vector.h"
#include "../algebraic/LaplacianMatrix.h"
#include "../properties/ConnectedComponents.h"
#include "../auxiliary/Timer.h"

#include "Config.h"
#include "SpanningTree.h"
#include "CycleDistribution.h"
#include "Flow.h"

namespace NetworKit {
namespace SDD {

/** @addtogroup sdd
 *  @{ */

/**
 * @brief Solver for a linear system \f$Lx = b\f$ where
 * \f$L\f$ is the Laplacian matrix of a graph.
 */
/**
 * Solves the Laplacian system \f$L(G)x = b\f$.
 *
 * @copydetails solveSDD
 *
 * @a bug If you use LogFlow, the LCA data structure is computed twice.
 * Can this be fixed without undue coupling?
 */
template<typename TCycle = UniformCycleDistribution, typename TFlow = TrivialFlow, SpanningTreeAlgo STAlgo = minDistanceST>
Vector solveLaplacian(const Graph& G, const Vector& b, SolverStatus& status /* inout */);

/** @} */

namespace {
  // Primal energy of the given flow is sum of f_uv^2*r_uv
  template<typename TFlow>
  edgeweight computeEnergy(const RootedTree& T, const TFlow& flow,
                           const std::vector<edgeweight>& off_tree_flow, const std::vector<edgeweight>& off_tree_resistance) {
    edgeweight out = 0.0;
    for (index i = 0; i < off_tree_flow.size(); ++i) {
      out += std::pow(off_tree_flow[i], 2) * off_tree_resistance[i];
    }
    T.forWeightedEdges([&] (node u, node v, edgeweight ew) {
      out += std::pow(flow.getPotential(u, v), 2) * ew;
    });
    return out;
  }

  /* Complete solver for one connected component. */
  template <typename TCycle, typename TFlow, SpanningTreeAlgo STAlgo>
  Vector perComponent(Graph& G, const Vector& b, SolverStatus& status, SeedType cur_seed) {
    assert(G.numberOfNodes() == b.getDimension());
    assert(abs(vectorSum(b)) < EPSILON);
    count m = G.numberOfEdges(), n = G.numberOfNodes();
    Aux::Timer timer;
    using namespace std;

    /** BEGIN: init */
    timer.start();

    /* Choose root node randomly. */
    RandomEngine engine(cur_seed);
    auto rand_node = std::uniform_int_distribution<node>(0, G.numberOfNodes() - 1);
    node root = rand_node(engine);

    /* Compute spanning tree */
    auto T = STAlgo(G, root);

    /* Initialise off-tree-edges. */
    vector<Edge> off_tree_edges;
    vector<edgeweight> off_tree_resistances;
    G.forEdges([&] (node u, node v, edgeweight ew) {
      if (!T.hasEdge(u, v)) {
        off_tree_edges.emplace_back(u, v);
        off_tree_resistances.emplace_back(1. / ew);
      }
    });
    vector<edgeweight> off_tree_flow(off_tree_edges.size());

    /* Compute stretches and resistances of off-tree-edges */
    vector<edgeweight> off_tree_cycle_resistances = computeResistances(T, off_tree_edges);
    vector<edgeweight> stretches(off_tree_edges.size());
    edgeweight summed_stretch = n - 1;
    for (index i = 0; i < off_tree_edges.size(); ++i) {
      stretches[i] = off_tree_cycle_resistances[i] * off_tree_resistances[i];
      summed_stretch += stretches[i];
      off_tree_cycle_resistances[i] += off_tree_resistances[i];
    }

    /* Initialize Laplacian, cycle and flow data structures. */
    TFlow flow(T, b);
    if (off_tree_edges.size() == 0) { // If G is a tree, we are done
      return flow.primalSolution();
    }
    TCycle cycles(G, T, off_tree_edges, stretches, status.seed);

    /* Set flow to the flow given by the initial vector
     * plus an adjustment to get a feasible flow, if necessary */
    if (status.debug_enable_initial) {
      // Compute necessary remainder in order to fix the vector
      Vector remainder = b;
      T.forEdges([&] (node u, node v) {
        remainder[u] += status.debug_initial[u] - status.debug_initial[v];
        remainder[v] -= status.debug_initial[u] - status.debug_initial[v];
      });
      for (index i = 0; i < off_tree_edges.size(); ++i) {
        Edge e = off_tree_edges[i];
        remainder[e.u] += status.debug_initial[e.u] - status.debug_initial[e.v];
        remainder[e.v] -= status.debug_initial[e.u] - status.debug_initial[e.v];
      }

      // Set initial flow
      flow = TFlow(T, remainder);
      T.forEdges([&] (node u, node v) {
        flow.addFlow(u, v, status.debug_initial[v] - status.debug_initial[u]);
      });
      for (index i = 0; i < off_tree_edges.size(); ++i) {
        Edge e = off_tree_edges[i];
        off_tree_flow[i] = status.debug_initial[e.v] - status.debug_initial[e.u];
      }
    }

    timer.stop();
    status.time_init += timer.elapsedMilliseconds();
    /** END: init */

    /* Compute tree condition number and iteration number described in paper. */
    LaplacianMatrix L(G);
    double condition = summed_stretch + m - 2*n + 2;
    status.stretch += summed_stretch / m;
    if (condition > 0 && summed_stretch > 0) {
      status.computed_niters += int(floor(condition * log(summed_stretch*condition/(status.desired_residual*b.length()))));
    }

    /** BEGIN: main */
    timer.start();

    /* Determine tree scaling factors */
    vector<double> kappas;
    double prod_kappas = 1.;
    if (status.precondition) {
      double kappa = log(static_cast<double>(n));
      while (kappa > 2) {
        prod_kappas *= kappa;
        kappas.emplace_back(kappa);
        kappa = log(kappa);
      }
      flow.scaleResistances(1./prod_kappas);
    }
    kappas.emplace_back(1.);
    count ncondition = kappas.size();

    /* For all tree scalings */
    auto cur_residual = numeric_limits<double>::infinity();
    int residual_update_freq = status.residual_update_freq == 0 ? n : status.residual_update_freq;
    for (index i = 0; i < ncondition; ++i) {
      if (status.precondition) {
        flow.scaleResistances(kappas[i]);
      }

      /* Main loop of the algorithm: Repairs a cycle each iteration. */
      while (cur_residual > status.desired_residual/ncondition && status.niters <= status.max_iters) {
        // Trace computed primal vectors if enable
        if (status.debug_enable_vectors) {
          status.debug_vectors.emplace_back(flow.primalSolution());
        }

        /* Get random cycle represented by e. */
        index e_idx = cycles.randomCycle();
        Edge e = off_tree_edges[e_idx];
        double& edge_flow = off_tree_flow[e_idx];

        /* Repair cycle e. */
        //TRACE("Potential of edge (", e.u, ", ", e.v, "): ", flow.getPotential(e.u), " ", flow.getPotential(e.v));
        double delta_flow = (edge_flow * off_tree_resistances[e_idx] - flow.getPotential(e.u, e.v)) / off_tree_cycle_resistances[e_idx];
        flow.addFlow(e.u, e.v, delta_flow);
        edge_flow += -delta_flow;

        /* DEBUG: length of cycle vs energy */
        if (status.debug_enable_stretch_vs_energy) {
          double a_opt = (edge_flow * off_tree_resistances[e_idx] - flow.getPotential(e.u, e.v));
          status.debug_stretch.emplace_back(off_tree_cycle_resistances[e_idx] / off_tree_resistances[e_idx]);
          status.debug_stretch_energy.emplace_back(a_opt * a_opt / off_tree_cycle_resistances[e_idx]);
          status.debug_stretch_time.emplace_back(flow.getTime(e.u) + flow.getTime(e.v));
        }

        /* Now update the residual (heuristically determine how often to do so) */
        if (status.niters % residual_update_freq == 0) {
          Vector x = flow.primalSolution();
          cur_residual = residual(L, x, b);

          if (status.debug_trace) {
            status.debug_residuals.emplace_back(cur_residual);
            status.debug_energy.emplace_back(computeEnergy(T, flow, off_tree_flow, off_tree_resistances));
            status.debug_iters.emplace_back(status.niters);
            status.debug_times.emplace_back(timer.elapsedMilliseconds());
            if (timer.elapsedMilliseconds() > status.max_time) {
              break;
            }
          }
#if (PRINT_RESIDUAL)
          std::cout << std::fixed << std::setprecision(10);
          Aux::printToStreamF(std::cout, "  Residual: %s, Energy: %s, Iters: %s\n",
                              cur_residual, computeEnergy(T, flow, off_tree_flow, off_tree_resistances), status.niters);
          std::cout.flush();
#endif
        }

        ++status.niters;
      }
    }

    timer.stop();
    status.time_main += timer.elapsedMilliseconds();
    /** END: main */

    status.actual_residual += cur_residual;
    return flow.primalSolution();
  }
}

template <typename TCycle, typename TFlow, SpanningTreeAlgo STAlgo>
Vector solveLaplacian(const Graph &G, const Vector& b, SolverStatus& status /* inout */) {
  using namespace std;
  if (G.numberOfNodes() == 0) {
    throw std::invalid_argument("empty graph disallowed");
  }
  if (G.numberOfNodes() != b.getDimension()) {
    throw std::invalid_argument("Dimension of b should be the order of G");
  }
  assert(status.desired_residual >= 0.);

  /* Get the components */
  TRACE("decomposing G into components");
  ConnectedComponents con(G);
  con.run();
  TRACE("Number of components:", con.numberOfComponents());
  unordered_map<index, vector<node>> comp_to_nodes;
  G.forNodes([&] (node v) {
    comp_to_nodes[con.componentOfNode(v)].emplace_back(v);
  });

  /* Build seed sequence from original seed */
  RandomEngine seeds(status.seed);

  /* For each component: run main algorithm */
  Vector x(G.numberOfNodes(), 0.0);
  double residual = status.desired_residual;
  status.niters = 0;
  status.time_init = 0;
  status.time_main = 0;
  status.actual_residual = 0.0;
  status.stretch = 0.0;
  status.computed_niters = 0;
  status.debug_residuals.clear();
  status.debug_times.clear();

//  cout << comp_to_nodes.size() << endl;
  for (const auto& comp: comp_to_nodes) {
    const auto& idx = comp.second;
    Graph Gcomp = inducedSubgraph(G, idx);
    Vector bcomp = inducedVector(b, idx);

    double sum = vectorSum(bcomp);
//    cout << sum << endl;
    if (abs(sum) > EPSILON) {
      for (index i = 0; i < bcomp.getDimension(); ++i) {
        bcomp[i] -= sum / bcomp.getDimension();
      }
//      WARN("Sum of entries of b should be 0 in every component");
    }

    status.desired_residual = residual * Gcomp.numberOfNodes() / G.numberOfNodes();
    Vector xcomp = perComponent<TCycle, TFlow, STAlgo>(Gcomp, bcomp, status, seeds());

    /* Translate solution back */
    for (index i = 0; i < xcomp.getDimension(); ++i) {
      x[idx[i]] = xcomp[i];
    }
  }

  status.converged = status.niters <= status.max_iters;
  status.desired_residual = residual;

  return x;
}

}
}

#endif
