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

Vector GaussSeidelRelaxation::relax(const CSRMatrix &A, const Vector &b, const Vector &initialGuess, const count maxIterations) const {
	count iterations = 0;
	Vector x_old = initialGuess;
	Vector x_new = initialGuess;
	if (maxIterations == 0) return initialGuess;

	count dimension = A.numberOfColumns();
	Vector diagonal = A.diagonal();

	do {
		x_old = x_new;

		for (index i = 0; i < dimension; ++i) {
			double sigma = 0.0;
			A.forNonZeroElementsInRow(i, [&](index column, double value) {
				if (column != i) {
					sigma += value * x_new[column];
				}
			});

			x_new[i] = (b[i] - sigma) / diagonal[i];
		}

		iterations++;
	} while (iterations < maxIterations);

	return x_new;
}

Vector GaussSeidelRelaxation::relax(const CSRMatrix &A, const Vector &b, const count maxIterations) const {
	Vector x(b.getDimension());
	return relax(A, b, x, maxIterations);
}


} /* namespace NetworKit */
