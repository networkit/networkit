/*
 * LaplacianMatrix.cpp
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "LaplacianMatrix.h"

namespace NetworKit {

LaplacianMatrix::LaplacianMatrix(const Graph &graph) : Matrix(graph.upperNodeIdBound()) {
	std::vector<std::pair<index, index>> positions;
	std::vector<double> values;

	graph.forNodes([&](const index i){
		double weightedDegree = graph.weightedDegree(i);

		double selfLoopWeight = 0.0;
		graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
			if (j == i) {
				selfLoopWeight = weight;
			} else {
				positions.push_back(std::make_pair(i,j));
				values.push_back(-weight);
			}
		});

		positions.push_back(std::make_pair(i,i));
		values.push_back(weightedDegree - selfLoopWeight); // degree matrix
	});

	graph.forNodes([&](const index i){
		double weightedDegree = graph.weightedDegree(i);

		graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
			setValue(i, j, -weight);
		});

		setValue(i, i, weightedDegree - graph.weight(i, i)); // degree matrix
	});
}


} /* namespace NetworKit */
