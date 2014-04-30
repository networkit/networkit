/*
 * AdjacencyMatrix.cpp
 *
 *  Created on: 28.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "AdjacencyMatrix.h"

namespace NetworKit {

AdjacencyMatrix::AdjacencyMatrix(const Graph &graph) : Matrix(graph.numberOfNodes()) {
	graph.forEdges([&](const node &i, const node &j) {
		setValue(i, j, 1.0);
	});
}

AdjacencyMatrix::~AdjacencyMatrix() {
}



} /* namespace NetworKit */

