/*
 * LinearSolver.h
 *
 *  Created on: 30.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_

#include "../algebraic/Matrix.h"
#include "../algebraic/LaplacianMatrix.h"
#include "../algebraic/Vector.h"
#include "../graph/Graph.h"

namespace NetworKit {

class LinearSolver {
protected:
	double tolerance;

public:
	LinearSolver(const double tolerance);
	virtual ~LinearSolver();

	virtual void setup(const Matrix &matrix) = 0;
	void setup(const Graph &graph);


	virtual bool solve(const Vector &rhs, Vector &result, count maxIterations = 20) = 0;

	bool isConverged(const Vector &residual) const;
};

} /* namespace NetworKit */

#endif /* LINEARSOLVER_H_ */
