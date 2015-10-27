///*
// * SparsifySolver.cpp
// *
// *  Created on: 23.06.2014
// *      Author: dhoske
// */
//
//#ifndef NOEIGEN
//
//#include "Config.h"
//#include "SparsifySolver.h"
//#include "SpanningTree.h"
//#include "../graph/GraphBuilder.h"
//
//using namespace std;
//
//namespace NetworKit {
//namespace SDD {
//
//SparsifySolver::SparsifySolver(const EigenMatrix& M, count depth)
//    : depth(depth), Dinvs(depth + 1), I_ADinvs(depth + 1), I_DinvAs(depth + 1) {
//  //assert(isSDD(Min));
//  //assert(isConnected(Min));
//
//  count n = M.rows();
//  EigenMatrix I(n, n);
//  I.setIdentity();
//
//  // Decompose M into M = D - A
//  EigenMatrix Dcur(n, n);
//  // Does not have predefined constructor for extracting a diagonal matrix?
//  for (index j = 0; j < n; ++j) {
//    Dcur.insert(j, j) = M.coeff(j, j);
//  }
//  EigenMatrix Acur = Dcur - M;
//
//  // Build matrices whose product approximates M^-1
//  for (index i = 0; i <= depth; ++i) {
//    // Store D^-1, I - AD^-1 and I - D^-1A
//    Dinvs[i] = Dcur.cwiseInverse();
//    EigenMatrix ADinv = Acur*Dinvs[i];
//    I_ADinvs[i] = I + ADinv;
//    I_DinvAs[i] = I_ADinvs[i].transpose();
//
//    Acur = ADinv * Acur;
//    /** @todo Sparsify Mcur = Dcur - Acur */
//    // Avoid the multiplication A*A that generates a dense matrix
//    //Acur = parallelSparsify(Acur, 1./depth, depth);
//  }
//}
//
//EigenVector SparsifySolver::run(const EigenVector& b) const {
//  // Compute result b forwards
//  std::vector<EigenVector> bs;
//  bs.emplace_back(b);
//  for (index i = 1; i <= depth; ++i) {
//    bs.emplace_back(I_ADinvs[i - 1] * bs[i - 1]);
//  }
//
//  // Compute solution x backwards
//  EigenVector x = Dinvs[depth] * bs[depth];
//  for (index i = depth - 1; i != index(-1); --i) {
//    x = 0.5 * (Dinvs[i] * bs[i] + I_DinvAs[i] * x);
//  }
//  return x;
//}
//
//namespace {
//  // A single round of the parallel sparsification
//  Graph parallelSparsifyRound(const Graph& G, double epsilon, default_random_engine& engine) {
//    auto unit = uniform_real_distribution<edgeweight>(0.0, 1.0);
//    //double tmp = log(G.numberOfNodes()) / epsilon;
//    //count t = static_cast<count>(floor(24*tmp*tmp));
//    count t = 1;
//
//    // Add spanning trees.
//    /** @todo: parallelise spanner construction */
//    Graph Gcopy = G;
//    GraphBuilder Gout(G.numberOfNodes(), true);
//    for (index i = 0; i < t && isConnected(Gcopy); ++i) {
//      auto T = minDistanceST(Gcopy, 0);
//      T.forWeightedEdges([&] (node u, node v, edgeweight w) {
//        Gcopy.removeEdge(u, v);
//        Gout.addEdge(u, v, w);
//      });
//    }
//
//    // Deal with remaining edges in parallel
//    G.parallelForEdges([&] (node u, node v, edgeweight w) {
//      if (unit(engine) < 1./4.) {
//        Gout.addEdge(u, v, w);
//      }
//    });
//
//    return Gout.toGraph();
//  }
//
//}
//
//EigenMatrix parallelSparsify(const EigenMatrix& M, double epsilon, double rho, SeedType seed) {
//  assert(epsilon > 0 && rho > 0);
//  default_random_engine engine(seed);
//
//  // Turn into a graph. Caveat: its pretty silly to convert between matrices and graphs
//  count n = M.rows();
//  Graph G(n, true, true);
//  for (int k = 0; k < M.outerSize(); ++k) {
//    for (EigenMatrix::InnerIterator it(M, k); it; ++it) {
//      G.addEdge(it.row(), it.col(), it.value());
//    }
//  }
//
//  // Sparsify
//  const count rounds = static_cast<count>(ceil(log(rho)));
//  for (index i = 0; i < rounds; ++i) {
//    // Warning: G is changed in call to parallelSparsify
//    G = parallelSparsifyRound(G, epsilon / rounds, engine);
//  }
//
//  // Turn into matrix again.
//  EigenMatrix Mout(n, n);
//  G.forEdges([&] (node u, node v, edgeweight w) {
//    Mout.coeffRef(u, v) = w;
//  });
//  return Mout;
//}
//
//}
//}
//
//#endif
