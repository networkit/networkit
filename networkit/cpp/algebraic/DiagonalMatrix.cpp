/*
 * DiagonalMatrix.cpp
 *
 *  Created on: 13.11.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "DiagonalMatrix.h"

namespace NetworKit {

DiagonalMatrix::DiagonalMatrix(const count dimension, const std::vector<double> &values) : Matrix(dimension) {
	for (index i = 0; i < dimension; ++i) {
		setValue(i, i, values[i]);
	}
}

} /* namespace NetworKit */
