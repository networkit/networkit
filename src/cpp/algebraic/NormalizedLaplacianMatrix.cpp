/*
 * NormalizedLaplacianMatrix.cpp
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "NormalizedLaplacianMatrix.h"

NormalizedLaplacianMatrix::NormalizedLaplacianMatrix() : Matrix() {
}

NormalizedLaplacianMatrix::NormalizedLaplacianMatrix(const NetworKit::Graph &graph) : Matrix(graph.numberOfNodes()) {
	auto rowIterator = [&](const uint64_t &i) {
			auto neighborIterator = [&](const uint64_t &j) {
				if (i != j) {
					uint64_t iDegree = graph.degree(i);
					uint64_t jDegree = graph.degree(j);
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
