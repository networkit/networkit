/*
 * IncidenceMatrix.cpp
 *
 *  Created on: 21.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "IncidenceMatrix.h"

namespace NetworKit {

IncidenceMatrix::IncidenceMatrix(const Graph &graph) : Matrix(graph.upperNodeIdBound(), graph.numberOfEdges()) {
	if (graph.isDirected()) {
		index edgeId = 0;
		graph.forEdges([&](node u, node v, edgeweight weight) {
			if (u != v) {
				edgeweight w = sqrt(weight);
				setValue(u, edgeId, w);
				setValue(v, edgeId, -w);
			}
			edgeId++;
		});
	} else {
		index edgeId = 0;
		graph.forEdges([&](node u, node v, edgeweight weight){
			if (u != v) {
				edgeweight w = sqrt(weight);
				if (u < v) { // orientation: small node number -> great node number
					setValue(u, edgeId, w);
					setValue(v, edgeId, -w);
				} else {
					setValue(u, edgeId, -w);
					setValue(v, edgeId, w);
				}
			}
			edgeId++;
		});
	}
}


} /* namespace NetworKit */





