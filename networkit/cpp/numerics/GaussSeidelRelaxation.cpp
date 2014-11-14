/*
 * GaussSeidelRelaxation.cpp
 *
 *  Created on: 27.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "GaussSeidelRelaxation.h"

namespace NetworKit {

GaussSeidelRelaxation::GaussSeidelRelaxation(double tolerance) : tolerance(tolerance) {
}

Vector GaussSeidelRelaxation::relax(const Matrix &A, const Vector &b, const Vector &initialGuess, const count maxIterations) const {
	count iterations = 0;
	Vector x_old = initialGuess;
	Vector x_new = initialGuess;
	do {
		x_old = x_new;

		for (index i = 0; i < x_new.getDimension(); ++i) {
			double sigma = 0.0;
			A.forNonZeroElementsInRow(i, [&](index row, index column, double value){
				if (column != row) {
					sigma += value * x_new[column];
				}
			});

			x_new[i] = (b[i] - sigma) / A(i,i);
		}


		iterations++;
	} while ((x_new - x_old).length() >= tolerance && iterations < maxIterations);

	return x_new;
}


} /* namespace NetworKit */
