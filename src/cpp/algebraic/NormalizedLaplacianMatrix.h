/*
 * NormalizedLaplacianMatrix.h
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NORMALIZEDLAPLACIANMATRIX_H_
#define NORMALIZEDLAPLACIANMATRIX_H_

#include "Matrix.h"
#include "../graph/Graph.h"
#include <math.h>

class NormalizedLaplacianMatrix : public Matrix {
public:
	NormalizedLaplacianMatrix();

	/**
	 * Constructs a normalized Laplacian matrix for the given @a graph.
	 */
	NormalizedLaplacianMatrix(const NetworKit::Graph &graph);
	virtual ~NormalizedLaplacianMatrix();
};

#endif /* NORMALIZEDLAPLACIANMATRIX_H_ */
