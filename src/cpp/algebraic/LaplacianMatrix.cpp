/*
 * LaplacianMatrix.cpp
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "LaplacianMatrix.h"

namespace NetworKit {

LaplacianMatrix::LaplacianMatrix() : Matrix() {
}

LaplacianMatrix::LaplacianMatrix(const Graph &graph) : Matrix(graph.numberOfNodes()) {
	auto rowIterator = [&](const index &i) {
		setValue(i, i, graph.degree(i));
		auto neighborIterator = [&](const index &j) {
			if (i != j) {
				setValue(i, j, -1.0);
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



} /* namespace NetworKit */
