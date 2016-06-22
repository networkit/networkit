/*
 * AdjacencyMatrix.cpp
 *
 *  Created on: 28.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "AdjacencyMatrix.h"

namespace NetworKit {

AdjacencyMatrix::AdjacencyMatrix(const Graph &graph) : Matrix(graph.upperNodeIdBound()) {
	graph.forEdges([&](const node &i, const node &j, double edgeWeight) {
		setValue(i, j, edgeWeight);
		if (!graph.isDirected()) {
			setValue(j, i, edgeWeight);
		}
	});
}


} /* namespace NetworKit */

