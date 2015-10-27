/*
 * DiagonalMatrix.cpp
 *
 *  Created on: 13.11.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "DiagonalMatrix.h"

namespace NetworKit {

DiagonalMatrix::DiagonalMatrix(const count dimension, const std::vector<double> &values) : Matrix(dimension) {
	if (values.size() != dimension) {
		throw std::runtime_error("DiagonalMatrix::DiagonalMatrix(count dimension, std::vector<double> values): dimension of values does not match the specified dimension");
	}

	for (index i = 0; i < dimension; ++i) {
		setValue(i, i, values[i]);
	}
}

} /* namespace NetworKit */
