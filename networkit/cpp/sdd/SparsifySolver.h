/*
 * SparsifySolver.h
 *
 *  Created on: 23.06.2014
 *      Author: dhoske
 */

//#ifndef NOEIGEN
//#ifndef SPARSIFY_SOLVER_H
//#define SPARSIFY_SOLVER_H
//
//#include <vector>
//
//#include "Config.h"
//#include "../algebraic/EigenLib.h"
//
//namespace NetworKit {
//namespace SDD {
//
///**
// * SDD-solver based on sparsification. (See: http://arxiv.org/pdf/1311.3286.pdf).
// * INCOMPLETE
// */
//class SparsifySolver {
//public:
//  SparsifySolver(const EigenMatrix& A, count depth);
//  EigenVector run(const EigenVector& b) const;
//
//  /* Neither copyable nor moveable */
//  SparsifySolver(const SparsifySolver&) = delete;
//  SparsifySolver(SparsifySolver&&) = delete;
//  SparsifySolver& operator=(const SparsifySolver&) = delete;
//  SparsifySolver& operator=(SparsifySolver&&) = delete;
//
//private:
//  count depth;
//  std::vector<EigenMatrix> Dinvs;
//  std::vector<EigenMatrix> I_ADinvs;
//  std::vector<EigenMatrix> I_DinvAs;
//};
//
///**
// * Parallel sparsification routine that returns a matrix \f$\tilde{M}\f$
// * fulfilling \f$(1-\epsilon)M \leq \tilde{M} \leq (1+\epsilon)M \f$ with high probability.
// * The expected number of its edges is \f$O(n\log(n)^3\log(rho)^3/\epsilon^2 + m/\rho)\f$.
// */
//EigenMatrix parallelSparsify(const EigenMatrix& M, double epsilon, double rho, SeedType seed = DEFAULT_SEED);
//
//}
//}
//
//#endif
//#endif
