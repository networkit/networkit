/*
 * LAMGSolver.h
 *
 *  Created on: 25.11.2014
 *      Author: Michael
 */

#ifndef LAMGSOLVER_H_
#define LAMGSOLVER_H_

#include "LinearSolver.h"
#include "MultigridHierarchyBuilder.h"
#include "MultigridHierarchy.h"
#include "Smoother.h"
#include "CG.h"
#include "Preconditioner.h"
#include "LAMGHierarchy.h"


namespace NetworKit {

class LAMGSolver : public LinearSolver {
private:
	const MultigridHierarchyBuilder<LAMGHierarchy> &builder;
	LAMGHierarchy hierarchy;
	const Smoother &smoother;

	Vector multigridCycle(const index level, const Vector &currentApproximation, const Vector &rhs);
//	Vector cycle(const Vector &xFine, const Vector &bFine);

public:
	LAMGSolver(const double tolerance, const MultigridHierarchyBuilder<LAMGHierarchy> &builder, const Smoother &smoother);

	using LinearSolver::setup;
	void setup(const Matrix &matrix);

	bool solve(const Vector &rhs, Vector &result, count maxIterations = 20);
	bool solve(const std::string &graph_name, const Vector &rhs, Vector &result, count maxIterations = 20);
};

} /* namespace NetworKit */

#endif /* LAMGSOLVER_H_ */
