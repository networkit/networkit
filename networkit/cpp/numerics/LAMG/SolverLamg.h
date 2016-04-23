/*
 * SolverLamg.h
 *
 *  Created on: 12.01.2015
 *      Author: Michael
 */

#ifndef SOLVERLAMG_H_
#define SOLVERLAMG_H_

#include "LevelHierarchy.h"
#include "../Smoother.h"
#include "../../algebraic/DenseMatrix.h"

namespace NetworKit {

/**
 * Status parameters of the solver.
 */
struct LAMGSolverStatus {
	// in
	count maxIters = std::numeric_limits<count>::max(); // maximum number of iterations
	count maxConvergenceTime = std::numeric_limits<count>::max(); // maximum time in milliseconds spent to solve the system
	double desiredResidualReduction = 1e-8; // desired reduction of the initial residual (finalResidual <= desiredResReduction * initialResidual)
	count numPreSmoothIters = 1; // number of pre smoothing iterations
	count numPostSmoothIters = 2; // number of post smoothing iterations

	// out
	count numIters; // number of iterations needed during solve phase
	double residual; // absolute final residual
	bool converged; // flag of conversion status
	std::vector<double> residualHistory; // history of absolute residuals
};

/**
 * @ingroup numerics
 * Implements the solve phase of LAMG (Lean Algebraic Multigrid by Livne et al.).
 */
class SolverLamg {
private:
	LevelHierarchy &hierarchy;
	const Smoother &smoother;
#ifndef NDEBUG
	static count minResTime;
	static count interpolationTime;
	static count restrictionTime;
	static count coarsestSolve;
#endif

	// data structures for iterate recombination
	std::vector<std::vector<Vector>> history;
	std::vector<std::vector<Vector>> rHistory;
	std::vector<index> latestIterate;
	std::vector<count> numActiveIterates;

	// bStages for Elimination Levels
	std::vector<std::vector<Vector>> bStages;

	void solveCycle(Vector &x, const Vector &b, int finest, LAMGSolverStatus &status);
	void cycle(Vector &x, const Vector &b, int finest, int coarsest, std::vector<count> &numVisits, std::vector<Vector> &X, std::vector<Vector> &B, const LAMGSolverStatus &status);
	void multigridCycle(index level, Vector &xf, const Vector &bf);
	void saveIterate(index level, const Vector &x, const Vector &r);
	void clearHistory(index level);
	void minRes(index level, Vector &x, const Vector &r);

public:
	/**
	 * Constructs a new solver instance for the specified @a hierarchy. The @a smoother will be used for relaxing and
	 * solving the coarser solutions.
	 * @param hierarchy Reference to the LevelHierarchy constructed by MultiLevelSetup.
	 * @param smoother Reference to a smoother.
	 */
	SolverLamg(LevelHierarchy &hierarchy, const Smoother &smoother);

	SolverLamg (const SolverLamg &other) = default;

	SolverLamg (SolverLamg &&other) = default;

	virtual ~SolverLamg() = default;

	SolverLamg& operator=(SolverLamg &&other) = default;

	SolverLamg& operator=(const SolverLamg &other) = default;

	/**
	 * Solves the system A*x = b for the given initial @a x and right-hand side @a b. More parameters can be specified
	 * in @a status and additional output is also stored in @a status. After the solver finished, the approximate
	 * solution is stored in @a x.
	 * @param x[out] Reference to the initial guess to the solution and the approximation after the solver finished.
	 * @param b The right-hand side vector.
	 * @param status Reference to an LAMGSolverStatus.
	 */
	void solve(Vector &x, const Vector &b, LAMGSolverStatus &status);
};

} /* namespace NetworKit */

#endif /* SOLVERLAMG_H_ */
