/*
 * LaplacianMatrix.h
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef LAPLACIANMATRIX_H_
#define LAPLACIANMATRIX_H_

#include "Matrix.h"
#include "../graph/Graph.h"
#include <cmath>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Laplacian matrix of a Graph.
 */
class LaplacianMatrix : public Matrix {
public:
	/**
	 * Constructs the LaplacianMatrix for the given @a graph.
	 */
	LaplacianMatrix(const Graph &graph);

};


} /* namespace NetworKit */

#endif /* LAPLACIANMATRIX_H_ */
