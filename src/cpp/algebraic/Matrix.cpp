/*
 * Matrix.cpp
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "Matrix.h"

Matrix::Matrix() : graph(0, true) {
}

Matrix::Matrix(const uint64_t &dimension) : graph(dimension, true) {
}

Matrix::Matrix(const uint64_t &dimension, const std::vector<std::pair<int, int> > &positions,
				const std::vector<double> &values) : graph(dimension, true) {
	assert(positions.size() == values.size());

	for (size_t k = 0; k < positions.size(); ++k) {
		assert(positions.at(k).first >= 0 && positions.at(k).second >= 0);

		std::pair<uint64_t, uint64_t> pos = positions.at(k);
		assert(pos.first < dimension);
		if (!graph.hasEdge(pos.first, pos.second)) { // symmetric matrix
			graph.addEdge(pos.first, pos.second, values.at(k));
		}
	}
}

Matrix::Matrix(const uint64_t &dimension, const std::vector<std::pair<NetworKit::node, NetworKit::node> > &positions, const std::vector<double> &values) : graph(dimension, true) {
	assert(positions.size() == values.size());

	for (size_t k = 0; k < positions.size(); ++k) {
		std::pair<NetworKit::node, NetworKit::node> pos = positions.at(k);
		assert(pos.first < static_cast<NetworKit::node>(dimension));
		if (!graph.hasEdge(pos.first,pos.second)) { // symmetric matrix
			graph.addEdge(pos.first, pos.second, values.at(k));
		}
	}
}

Matrix::Matrix(const std::vector<Vector> &columns) : graph(columns.size(), true) {
	// check dimension
	uint64_t dimension = columns.size();
#pragma omp parallel for
	for (size_t j = 0; j < dimension; ++j) {
		if (columns.at(j).getDimension() != dimension) {
			throw std::runtime_error("row dimensions of one or more columns do not match");
		}
	}

	for (uint64_t j = 0; j < dimension; ++j) {
		for (uint64_t i = 0; i < dimension; ++i) {
			double value = columns.at(j)(i);
			if (value != 0.0 && !graph.hasEdge(i, j)) { // do not store 0 values, symmetric matrix
				graph.addEdge(i, j, value);
			}
		}
	}
}

Matrix::Matrix(const Matrix &other) : graph(other.graph) {
}

Matrix::~Matrix() {
}

uint64_t Matrix::numberOfRows() const {
	return graph.numberOfNodes();
}

uint64_t Matrix::numberOfColumns() const {
	return numberOfRows();
}

double Matrix::operator()(const uint64_t &i, const uint64_t &j) const {
	if (i >= numberOfRows() || j >= numberOfColumns()) {
		throw std::out_of_range("at least one index out of range");
	}

	return graph.weight(i,j);
}

Vector Matrix::row(const uint64_t &i) const {
	std::vector<double> values(numberOfColumns(), 0.0);
	auto setElements = [&](NetworKit::node i, NetworKit::node j, double value) {
		values[j] = value;
	};

	graph.forWeightedEdgesOf(i, setElements);

	return Vector(values);
}

Vector Matrix::column(const uint64_t &j) const {
	std::vector<double> values(numberOfRows(), 0.0);
	auto setElements = [&](NetworKit::node i) {
		values[i] = graph.weight(i, j);
	};

	graph.parallelForNodes(setElements);

	return Vector(values);
}

Matrix Matrix::operator+(const Matrix &other) const {
	return Matrix(*this) += other;
}

Matrix& Matrix::operator+=(const Matrix &other) {
	if (numberOfRows() != other.numberOfRows() || numberOfColumns() != other.numberOfColumns()) {
		throw std::runtime_error("dimensions of matrices do not match");
	}

	auto add = [&](NetworKit::node i, NetworKit::node j, double value) {
		if (i >= j) { // symmetric matrix
			graph.increaseWeight(i, j, value);
		}
	};

	other.parallelForNonZeroElementsInRowOrder(add);

	return *this;
}

Matrix Matrix::operator-(const Matrix &other) const {
	return Matrix(*this) -= other;
}

Matrix& Matrix::operator-=(const Matrix &other) {
	if (numberOfRows() != other.numberOfRows() || numberOfColumns() != other.numberOfColumns()) {
		throw std::runtime_error("dimensions of matrices do not match");
	}

	auto subtract = [&](NetworKit::node i, NetworKit::node j, double value) {
		if (i >= j) { // symmetric matrix
			graph.increaseWeight(i, j, -value);
		}
	};
	other.parallelForNonZeroElementsInRowOrder(subtract);

	return *this;
}

Matrix Matrix::operator*(const double &scalar) const {
	return Matrix(*this) *= scalar;
}

Matrix& Matrix::operator*=(const double &scalar) {
	auto multiplyElement = [&](NetworKit::node i, NetworKit::node j, double value) {
		graph.setWeight(i, j, value * scalar);
	};
	graph.parallelForWeightedEdges(multiplyElement);

	return *this;
}

Vector Matrix::operator*(const Vector &vector) const {
	if (numberOfColumns() != vector.getDimension()) {
		throw std::runtime_error("dimensions of matrix and vector do not match");
	}

	Vector result(numberOfRows(), 0.0);

	auto multiplyElement = [&](NetworKit::node i, NetworKit::node j, double value) {
		result(i) += value * vector(j);
	};
	parallelForNonZeroElementsInRowOrder(multiplyElement);

	return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
	if (numberOfColumns() != other.numberOfRows()) {
		throw std::runtime_error("dimensions of matrices do not match");
	}

	std::vector<Vector> columns(other.numberOfColumns());
#pragma omp parallel for
	for (uint64_t j = 0; j < other.numberOfColumns(); ++j) {
		columns[j] = (*this * other.column(j));
	}

	return Matrix(columns);
}

