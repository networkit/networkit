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
#include <cmath>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Normalized laplacian matrix of a Graph.
 */
class NormalizedLaplacianMatrix : public Matrix {
public:
	/**
	 * Constructs the NormalizedLaplacianMatrix for the given @a graph.
	 */
	NormalizedLaplacianMatrix(const Graph &graph);
};

} /* namespace NetworKit */

#endif /* NORMALIZEDLAPLACIANMATRIX_H_ */
