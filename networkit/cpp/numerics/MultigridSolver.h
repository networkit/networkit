/*
 * MultigridSolver.h
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef MULTIGRIDSOLVER_H_
#define MULTIGRIDSOLVER_H_

#include "LinearSolver.h"
#include "MultigridHierarchy.h"
#include "MultigridHierarchyBuilder.h"
#include "Smoother.h"
#include "CG.h"
#include "Preconditioner.h"


namespace NetworKit {

class MultigridSolver : public LinearSolver {
private:
	count numPreSmooth;
	count numPostSmooth;

	const MultigridHierarchyBuilder &builder;
	MultigridHierarchy hierarchy;
	const Smoother &smoother;

	Vector restrictResidual(const Matrix &matrix, const Matrix &interpolationMatrix, const Vector &result, const Vector &rhs);
	void interpolateError(const Matrix &interpolationMatrix, const Vector &error, Vector &result);

	Vector vCycle(const Matrix &matrix, index level, const Vector &rhs, const Vector &result);
	Vector fCycle(const Matrix &matrix, index level, const Vector &rhs, const Vector &result);

public:
	MultigridSolver(const double tolerance, const count numPreSmooth, const count numPostSmooth,
										const MultigridHierarchyBuilder &builder, const Smoother &smoother);

	using LinearSolver::setup;
	void setup(const Matrix &matrix);

	bool solve(const Vector &rhs, Vector &result, count maxIterations = 20);
};

} /* namespace NetworKit */

#endif /* MULTIGRIDSOLVER_H_ */
