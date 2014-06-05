/*
 * IncidenceMatrix.cpp
 *
 *  Created on: 21.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "IncidenceMatrix.h"

namespace NetworKit {

IncidenceMatrix::IncidenceMatrix(Graph &graph) : graph(graph) {
}

double IncidenceMatrix::value(const node &node, const Edge &edge) const {
	if (edge.first != edge.second) {
		if (node == edge.first) {
			if (node > edge.second) {
				return -1.0;
			} else {
				return 1.0;
			}
		} else if (node == edge.second) {
			if (node < edge.first) {
				return 1.0;
			} else {
				return -1.0;
			}
		}
	}

	return 0.0;
}

count IncidenceMatrix::numberOfRows() const {
	return graph.numberOfNodes();
}

count IncidenceMatrix::numberOfColumns() const {
	return graph.numberOfEdges();
}

double IncidenceMatrix::operator()(const index &i, const index &j) const {
	return value(i, graph.edges().at(j));
}

Vector IncidenceMatrix::row(const index &i) const {
	Vector vector(numberOfColumns(), 0.0);
	std::vector<Edge> edges = graph.edges();

#pragma omp parallel for
	for (count k = 0; k < edges.size(); ++k) {
		vector[k] = value(i, edges.at(k));
	}

	return vector;
}

Vector IncidenceMatrix::column(const index &j) const {
	Vector vector(numberOfRows(), 0.0);
	Edge edge = graph.edges().at(j);

	if (edge.first != edge.second) {
		if (edge.first > edge.second) {
			vector[edge.first] = -1.0;
			vector[edge.second] = 1.0;
		} else {
			vector[edge.first] = 1.0;
			vector[edge.second] = -1.0;
		}
	}

	return vector;
}

Vector IncidenceMatrix::operator*(const Vector &vector) const {
	if (vector.isTransposed() || numberOfColumns() != vector.getDimension()) {
		throw std::runtime_error("operator *: nonconformant arguments");
	}

	Vector result(numberOfRows(), 0.0);
	std::vector<Edge> edges = graph.edges();

#pragma omp parallel for
	for (count i = 0; i < numberOfRows(); ++i) {
		for (count k = 0; k < edges.size(); ++k) {
			result[i] += value(i, edges.at(k)) * vector[k];
		}
	}

	return result;
}


} /* namespace NetworKit */





