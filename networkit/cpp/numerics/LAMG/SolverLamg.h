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

namespace NetworKit {

struct LAMGSolverStatus {
	count maxIters = std::numeric_limits<count>::max();
	count maxConvergenceTime = std::numeric_limits<count>::max();
	double desiredResidual = 1e-8;
	count numPreSmoothIters = 1;
	count numPostSmoothIters = 2;

	count numIters;
	double residual;
	bool converged;
	std::vector<double> residualHistory;
};

class SolverLamg {
private:
	LevelHierarchy &hierarchy;
	const Smoother &smoother;
#ifndef NPROFILE
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

	void solveCycle(Vector &x, const Vector &b, int finest, LAMGSolverStatus &status);
	void cycle(Vector &x, const Vector &b, int finest, int coarsest, std::vector<count> &numVisits, std::vector<Vector> &X, std::vector<Vector> &B, const LAMGSolverStatus &status);
	void multigridCycle(index level, Vector &xf, const Vector &bf);
	void saveIterate(index level, const Vector &x, const Vector &r);
	void clearHistory(index level);
	void minRes(index level, Vector &x, const Vector &r);

public:
	SolverLamg(LevelHierarchy &hierarchy, const Smoother &smoother);

	void solve(Vector &x, const Vector &b, LAMGSolverStatus &status);
};

} /* namespace NetworKit */

#endif /* SOLVERLAMG_H_ */
