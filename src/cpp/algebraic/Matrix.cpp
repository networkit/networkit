/*
 * Matrix.cpp
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "Matrix.h"

Matrix::Matrix() {
}

Matrix::Matrix(const int &dimension) : graph(dimension, true) {
}

Matrix::Matrix(const int &dimension, const std::vector<std::pair<NetworKit::node, NetworKit::node> > &positions, const std::vector<double> &values) : graph(dimension, true) {
	assert(positions.size() == values.size());

	for (size_t k = 0; k < positions.size(); ++k) {
		std::pair<NetworKit::node, NetworKit::node> pos = positions.at(k);
		assert(pos.first < static_cast<NetworKit::node>(dimension));
		graph.addEdge(pos.first, pos.second, values.at(k));
	}
}

Matrix::Matrix(const Matrix &other) : graph(other.graph) {
}

Matrix::~Matrix() {
}

int Matrix::numberOfRows() const {
	return graph.numberOfNodes();
}

int Matrix::numberOfColumns() const {
	return numberOfRows();
}

double Matrix::operator()(const uint64_t &i, const uint64_t &j) const {
	return graph.weight(i,j);
}

Matrix Matrix::operator+(const Matrix &other) const {
	return Matrix(*this) += other;
}

Matrix& Matrix::operator+=(const Matrix &other) {
	if (numberOfRows() != other.numberOfRows() || numberOfColumns() != other.numberOfColumns()) {
		throw std::runtime_error("dimensions of matrices do not match");
	}

	auto addRow = [&](NetworKit::node i) {
		auto addValues = [&](NetworKit::node i, NetworKit::node j, double value) {
			if (i >= j) {
				graph.increaseWeight(i, j, value);
			}
		};

		other.graph.forWeightedEdgesOf(i, addValues);
	};

	graph.balancedParallelForNodes(addRow);
	return *this;
}
