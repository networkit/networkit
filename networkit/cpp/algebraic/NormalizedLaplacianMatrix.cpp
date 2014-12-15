/*
 * NormalizedLaplacianMatrix.cpp
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "NormalizedLaplacianMatrix.h"

namespace NetworKit {

NormalizedLaplacianMatrix::NormalizedLaplacianMatrix(const Graph &graph) : Matrix(graph.upperNodeIdBound()) {
	graph.forNodes([&](const node i){
		double weightedDegree = graph.weightedDegree(i);

		graph.forNeighborsOf(i, [&](const node j, double weight){
			if (i != j) {
				double weightedNeighborDegree = graph.weightedDegree(j);
				setValue(i, j, -weight/sqrt(weightedDegree * weightedNeighborDegree));
			}
		});

		if (weightedDegree != 0.0) {
			if (graph.isWeighted()) {
				setValue(i, i, 1-(graph.weight(i, i)) / weightedDegree);
			} else {
				setValue(i, i, 1);
			}
		}
	});
}

} /* namespace NetworKit */
