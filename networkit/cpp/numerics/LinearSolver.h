/*
 * LinearSolver.h
 *
 *  Created on: 30.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_

#include "../algebraic/CSRMatrix.h"
#include "../algebraic/LaplacianMatrix.h"
#include "../algebraic/Vector.h"
#include "../graph/Graph.h"
#include <limits>

namespace NetworKit {

struct SolverStatus {
	count numIters; // number of iterations needed during solve phase
	double residual; // absolute final residual
	bool converged; // flag of conversion status
};

class LinearSolver {
protected:
	double tolerance;

public:
	LinearSolver(const double tolerance);
	virtual ~LinearSolver();

	virtual void setup(const CSRMatrix &matrix) = 0;
	virtual void setup(const Graph &graph);


	virtual SolverStatus solve(const Vector &rhs, Vector &result, count maxConvergenceTime = 5 * 60 * 1000, count maxIterations = std::numeric_limits<count>::max()) = 0;
};

} /* namespace NetworKit */

#endif /* LINEARSOLVER_H_ */
