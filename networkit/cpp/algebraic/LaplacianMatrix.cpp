/*
 * LaplacianMatrix.cpp
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "LaplacianMatrix.h"

namespace NetworKit {

LaplacianMatrix::LaplacianMatrix(const Graph &graph) : Matrix(graph.upperNodeIdBound()) {
	graph.forNodes([&](const index i){
		double weightedDegree = graph.weightedDegree(i);

		graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
			setValue(i, j, -weight);
		});

		setValue(i, i, weightedDegree - graph.weight(i, i)); // degree matrix
	});
}


} /* namespace NetworKit */
