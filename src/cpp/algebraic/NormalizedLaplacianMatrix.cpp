/*
 * NormalizedLaplacianMatrix.cpp
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "NormalizedLaplacianMatrix.h"

namespace NetworKit {

NormalizedLaplacianMatrix::NormalizedLaplacianMatrix() : Matrix() {
}

NormalizedLaplacianMatrix::NormalizedLaplacianMatrix(const Graph &graph) : Matrix(graph.numberOfNodes()) {
	auto rowIterator = [&](const index &i) {
		auto neighborIterator = [&](const index &j) {
			if (i != j) {
				count iDegree = graph.degree(i);
				count jDegree = graph.degree(j);
				double value = -1.0 / (sqrt(iDegree * jDegree));

				setValue(i, j, value);
			}
		};
		graph.forNeighborsOf(i, neighborIterator);

		if (graph.degree(i) != 0) {
			setValue(i, i, 1.0);
		}
	};

	graph.forNodes(rowIterator);
}

NormalizedLaplacianMatrix::~NormalizedLaplacianMatrix() {
}

} /* namespace NetworKit */
