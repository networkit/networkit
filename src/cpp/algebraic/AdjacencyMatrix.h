/*
 * AdjacencyMatrix.h
 *
 *  Created on: 28.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef ADJACENCYMATRIX_H_
#define ADJACENCYMATRIX_H_

#include "Matrix.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup algebraic
 * Adjacency matrix of a Graph.
 */
class AdjacencyMatrix : public Matrix {
public:
	/**
	 * Constructs the AdjacencyMatrix of @a graph.
	 */
	AdjacencyMatrix(const Graph &graph);
};

} /* namespace NetworKit */

#endif /* ADJACENCYMATRIX_H_ */
