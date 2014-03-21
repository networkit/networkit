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

class LaplacianMatrix : public Matrix {
public:
	LaplacianMatrix();

	/**
	 * Constructs a Laplacian matrix for the given @a graph.
	 */
	LaplacianMatrix(const NetworKit::Graph &graph);
	virtual ~LaplacianMatrix();
};

#endif /* LAPLACIANMATRIX_H_ */
