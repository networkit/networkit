/*
 * IncidenceMatrix.cpp
 *
 *  Created on: 21.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "IncidenceMatrix.h"

IncidenceMatrix::IncidenceMatrix(NetworKit::Graph &graph) : graph(graph) {
}

IncidenceMatrix::~IncidenceMatrix() {
}

double IncidenceMatrix::value(const uint64_t &node, const Edge &edge) const {
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

uint64_t IncidenceMatrix::numberOfRows() const {
	return graph.numberOfNodes();
}

uint64_t IncidenceMatrix::numberOfColumns() const {
	return graph.numberOfEdges();
}

double IncidenceMatrix::operator()(const uint64_t &i, const uint64_t &j) const {
	return value(i, graph.edges().at(j));
}

Vector IncidenceMatrix::row(const uint64_t &i) const {
	Vector vector(numberOfColumns(), 0.0);
	std::vector<Edge> edges = graph.edges();

#pragma omp parallel for
	for (uint64_t k = 0; k < edges.size(); ++k) {
		vector(k) = value(i, edges.at(k));
	}

	return vector;
}

Vector IncidenceMatrix::column(const uint64_t &j) const {
	Vector vector(numberOfRows(), 0.0);
	Edge edge = graph.edges().at(j);

	if (edge.first != edge.second) {
		if (edge.first > edge.second) {
			vector(edge.first) = -1.0;
			vector(edge.second) = 1.0;
		} else {
			vector(edge.first) = 1.0;
			vector(edge.second) = -1.0;
		}
	}

	return vector;
}

Vector IncidenceMatrix::operator*(const Vector &vector) const {
	if (vector.isTransposed() || numberOfColumns() != vector.getDimension()) {
		throw std::runtime_error("operator *: nonconformant arguments"); // TODO
	}

	Vector result(numberOfRows(), 0.0);
	std::vector<Edge> edges = graph.edges();

#pragma omp parallel for
	for (uint64_t i = 0; i < numberOfRows(); ++i) {
		for (uint64_t k = 0; k < edges.size(); ++k) {
			result(i) += value(i, edges.at(k)) * vector(k);
		}
	}

	return result;
}





