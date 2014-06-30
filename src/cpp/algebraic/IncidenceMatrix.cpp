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

double IncidenceMatrix::value(const node &nd, const Edge &edge) const {
	node u = edge.first;
	node v = edge.second;

	if (u != v) {
		double weight = sqrt(graph.weight(u,v));

		if (nd == u) {
			if (nd > v) {
				return -weight;
			} else {
				return weight;
			}
		} else if (nd == v) {
			if (nd < u) {
				return weight;
			} else {
				return -weight;
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

	node u = edge.first;
	node v = edge.second;

	if (u != v) {
		double weight = sqrt(graph.weight(u, v));
		if (u > v) {
			vector[u] = -weight;
			vector[v] = weight;
		} else {
			vector[u] = weight;
			vector[v] = -weight;
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





