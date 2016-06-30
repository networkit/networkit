/*
 * LinearSolver.cpp
 *
 *  Created on: 30.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "LinearSolver.h"

namespace NetworKit {

LinearSolver::LinearSolver(double tolerance) : tolerance(tolerance){
}

LinearSolver::~LinearSolver(){
}

void LinearSolver::setup(const Graph &graph) {
	setup(CSRMatrix::graphLaplacian(graph));
}

void LinearSolver::parallelSolve(const std::vector<Vector> &rhs, std::vector<Vector> &results, count maxConvergenceTime, count maxIterations) {
	for (index i = 0; i < rhs.size(); ++i) {
		solve(rhs[i], results[i], maxConvergenceTime, maxIterations);
	}
}

} /* namespace NetworKit */
