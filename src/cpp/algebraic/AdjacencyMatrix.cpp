/*
 * AdjacencyMatrix.cpp
 *
 *  Created on: 28.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "AdjacencyMatrix.h"

AdjacencyMatrix::AdjacencyMatrix() {
}

AdjacencyMatrix::AdjacencyMatrix(const NetworKit::Graph &graph) : Matrix(graph.numberOfNodes()) {
	auto edgeIterator = [&](const uint64_t &i, const uint64_t &j) {
		setValue(i, j, 1.0);
	};

	graph.forEdges(edgeIterator);
}

AdjacencyMatrix::~AdjacencyMatrix() {
}

