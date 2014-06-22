/*
 * LaplaceSolver.h
 *
 *  Created on: May 03, 2014
 *      Author: dhoske
 */

#ifndef LAPLACE_SOLVER_H
#define LAPLACE_SOLVER_H

#include <valgrind/callgrind.h>
#include <vector>
#include <unordered_map>

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

/**
 * @brief Solver for a linear system \f$Lx = b\f$ where
 * \f$L\f$ is the Laplace matrix of a graph.
 */
template<typename TCycle = UniformCycleDistribution, typename TFlow = TrivialFlow,
         SpanningTreeAlgo STAlgo = minDistanceST>
class LaplaceSolver {
public:
  LaplaceSolver() {
  }

  /**
   * Solves the system \f$L(G)x = b\f$ where \f$L(G)\f$ is the Laplace
   * matrix of @a G.
   */
  Vector run(const Graph& G, const Vector& b, SolverStatus& conf /* inout */);

  /* Deleted copiers */
  LaplaceSolver(const LaplaceSolver&) = delete;
  LaplaceSolver(LaplaceSolver&&) = delete;
  LaplaceSolver& operator=(const LaplaceSolver&) = delete;
  LaplaceSolver& operator=(LaplaceSolver&&) = delete;

private:
  /* Solves the Laplacian system for a connected graph. */
  Vector perComponent(const Graph& G, const Vector& b, SolverStatus& conf);
};

template <typename TCycle, typename TFlow, SpanningTreeAlgo STAlgo>
Vector LaplaceSolver<TCycle,TFlow,STAlgo>::run(const Graph &G, const Vector& b, SolverStatus& conf /* inout */) {
  using namespace std;

  /* Get the components */
  TRACE("decomposing G into components");
  ConnectedComponents con(G);
  con.run();
  TRACE("Number of components:", con.numberOfComponents());
  unordered_map<index, vector<node>> comp_to_nodes;
  for (node v = 0; v < G.numberOfNodes(); ++v) {
    comp_to_nodes[con.componentOfNode(v)].emplace_back(v);
  }

  /* For each component: run main algorithm */
  Vector x(G.numberOfNodes(), 0.0);
  conf.niters = 0;
  conf.time_init = 0;
  conf.time_main = 0;
  conf.actual_residual = 0.0;
  conf.desired_residual /= con.numberOfComponents();
  conf.stretch = 0.0;
  conf.computed_niters = 0;
  for (const auto& comp: comp_to_nodes) {
    const auto& idx = comp.second;
    Graph Gcomp = inducedSubgraph(G, idx);
    Vector bcomp = inducedVector(b, idx);
    Vector xcomp = perComponent(Gcomp, bcomp, conf);
    /* Translate solution back */
    for (index i = 0; i < xcomp.getDimension(); ++i) {
      x[idx[i]] = xcomp[i];
    }
  }
  conf.converged = conf.niters <= conf.max_iters;
  conf.desired_residual *= con.numberOfComponents();

  return x;
}

/* Complete solver for one connected component. */
template <typename TCycle, typename TFlow, SpanningTreeAlgo STAlgo>
Vector LaplaceSolver<TCycle,TFlow,STAlgo>::perComponent(const Graph& G, const Vector& b, SolverStatus& conf) {
  assert(G.numberOfNodes() == b.getDimension());
  assert(abs(vectorSum(b)) < EPSILON);
  Aux::Timer timer;
  using namespace std;

  timer.start();
  auto T = STAlgo(G, 0);
  //ERROR("Stretch: ", computeStretch(G, T));
  count n = G.numberOfNodes(), m = G.numberOfEdges();

  /* Flow on off-tree edges stored in a vector. */
  vector<Edge> off_tree_edges;
  vector<double> off_tree_resistance;
  G.forWeightedEdges([&] (node u, node v, edgeweight value) {
    if (!T.hasEdge(u, v)) {
      off_tree_resistance.emplace_back(1.0 / value);
      off_tree_edges.push_back({u, v});
    }
  });
  vector<double> off_tree_flow(off_tree_edges.size());

  /* Corresponding Laplacian is used only for computing the residual */
  LaplacianMatrix L(G);
  TFlow flow(T, b);
  TCycle cycles(G, T, off_tree_edges);
  auto cur_residual = numeric_limits<double>::infinity();
  timer.stop();
  conf.time_init += timer.elapsedMilliseconds();

  /* Number of iterations given in the simple solver in the paper */
  CALLGRIND_STOP_INSTRUMENTATION;
  double stretch = computeStretch(G, T);
  double condition = stretch + m - 2*n + 2;
  conf.stretch += stretch;
  if (condition > 0 && stretch > 0) {
    conf.computed_niters += int(floor(condition * log(stretch*condition/(conf.desired_residual*b.length()))));
  }
  CALLGRIND_START_INSTRUMENTATION;

  /* Main loop of the algorithm that repairs a cycle each iteration. */
  timer.start();
  int residual_update_freq = conf.residual_update_freq == 0 ? n*log(n) : conf.residual_update_freq;
  if (off_tree_edges.size() > 0) {
    while (cur_residual > conf.desired_residual && conf.niters <= conf.max_iters) {
      /* Get random cycle represented by e. */
      int e_idx = cycles.randomCycle();
      Edge e = off_tree_edges[e_idx];
      double& edge_flow = off_tree_flow[e_idx];
      double  edge_resistance = off_tree_resistance[e_idx];
      double  resistance =
          edge_resistance
        + flow.getResistance(e.u)
        + flow.getResistance(e.v);

      /* Repair cycle e. */
      //TRACE("Potential of edge (", e.u, ", ", e.v, "): ", flow.getPotential(e.u), " ", flow.getPotential(e.v));
      double delta_flow = (
          edge_flow * edge_resistance
        - flow.getPotential(e.u)
        + flow.getPotential(e.v)
      ) / resistance;

      flow.addFlow(e.u, delta_flow);
      flow.addFlow(e.v, -delta_flow);
      edge_flow += -delta_flow;

      /* Now update the residual (heuristically determine how often to do so) */
      if (conf.niters % residual_update_freq == 0) {
        Vector x = flow.primalSolution();
        cur_residual = residual(L, x, b);
        //DEBUG("Residual: ", cur_residual);
      }

      ++conf.niters;
    }
  }
  timer.stop();
  conf.time_main += timer.elapsedMilliseconds();
  conf.actual_residual += cur_residual;
  return flow.primalSolution();
}

}
}

#endif
