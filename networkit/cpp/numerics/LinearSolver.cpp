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
	setup(LaplacianMatrix(graph));
}

bool LinearSolver::isConverged(const Vector &residual) const {
	return residual.length() <= tolerance;
}

} /* namespace NetworKit */
