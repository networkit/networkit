/*
 * DiagonalMatrix.h
 *
 *  Created on: 13.11.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef DIAGONALMATRIX_H_
#define DIAGONALMATRIX_H_

#include "Matrix.h"

namespace NetworKit {

class DiagonalMatrix : public Matrix {
public:
	DiagonalMatrix(const count dimension, const std::vector<double> &values);
	inline DiagonalMatrix(const count dimension) : DiagonalMatrix(dimension, std::vector<double>(dimension, 1.0)) {}
};

} /* namespace NetworKit */

#endif /* DIAGONALMATRIX_H_ */
