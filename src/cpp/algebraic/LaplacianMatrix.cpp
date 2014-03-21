/*
 * LaplacianMatrix.cpp
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "LaplacianMatrix.h"

LaplacianMatrix::LaplacianMatrix() : Matrix() {
}

LaplacianMatrix::LaplacianMatrix(const NetworKit::Graph &graph) : Matrix(graph.numberOfNodes()) {
	auto rowIterator = [&](const uint64_t &i) {
		setValue(i, i, graph.degree(i));
		auto neighborIterator = [&](const uint64_t &j) {
			if (i != j) {
				setValue(i, j, -1);
			} else {
				setValue(i, i, (*this)(i, i) - 1.0);
			}
		};
		graph.forNeighborsOf(i, neighborIterator);
	};

	graph.forNodes(rowIterator);
}

LaplacianMatrix::~LaplacianMatrix() {
}

