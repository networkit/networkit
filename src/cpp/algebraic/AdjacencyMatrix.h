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

class AdjacencyMatrix : public Matrix {
public:
	AdjacencyMatrix();
	AdjacencyMatrix(const NetworKit::Graph &graph);
	virtual ~AdjacencyMatrix();
};

#endif /* ADJACENCYMATRIX_H_ */
