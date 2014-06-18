/*
 * Matrix.cpp
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "Matrix.h"

namespace NetworKit {

Matrix::Matrix() : graph(0, true) {
}

Matrix::Matrix(const uint64_t &dimension) : graph(dimension, true) {
}

Matrix::Matrix(const uint64_t &dimension, const std::vector<std::pair<int, int>> &positions,
				const std::vector<double> &values) : graph(dimension, true) {
	assert(positions.size() == values.size());

	for (size_t k = 0; k < positions.size(); ++k) {
		assert(positions[k].first >= 0 && positions[k].second >= 0);

		std::pair<uint64_t, uint64_t> pos = positions[k];
		assert(pos.first < dimension);
		if (!graph.hasEdge(pos.first, pos.second)) { // symmetric matrix
			graph.addEdge(pos.first, pos.second, values[k]);
		}
	}
}

Matrix::Matrix(const uint64_t &dimension, const std::vector<std::pair<node, node>> &positions, const std::vector<double> &values) : graph(dimension, true) {
	assert(positions.size() == values.size());

	for (size_t k = 0; k < positions.size(); ++k) {
		std::pair<node, node> pos = positions[k];
		assert(pos.first < dimension);
		if (!graph.hasEdge(pos.first,pos.second)) { // symmetric matrix
			graph.addEdge(pos.first, pos.second, values[k]);
		}
	}
}

Matrix::Matrix(const std::vector<Vector> &rows) : graph(rows.size(), true) {
	// check dimension
	uint64_t dimension = rows.size();
#pragma omp parallel for
	for (size_t i = 0; i < dimension; ++i) {
		if (rows[i].getDimension() != dimension) {
			throw std::runtime_error("row dimensions of one or more columns do not match");
		}
	}

	for (uint64_t i = 0; i < dimension; ++i) {
		for (uint64_t j = 0; j < dimension; ++j) {
			double value = rows[i][j];
			if (value != 0.0 && !graph.hasEdge(i, j)) { // do not store 0 values, symmetric matrix
				graph.addEdge(i, j, value);
			}
		}
	}
}

Matrix::Matrix(const Matrix &other) : graph(other.graph) {
}

uint64_t Matrix::numberOfRows() const {
	return graph.numberOfNodes();
}

uint64_t Matrix::numberOfColumns() const {
	return numberOfRows();
}

double Matrix::operator()(const index &i, const index &j) const {
	if (i >= numberOfRows() || j >= numberOfColumns()) {
		throw std::out_of_range("at least one index out of range");
	}

	return graph.weight(i,j);
}

void Matrix::setValue(const index &i, const index &j, const double &value) {
	graph.setWeight(i, j, value);
}

Vector Matrix::row(const index &i) const {
	Vector row(numberOfColumns(), 0.0, true);
	graph.forWeightedEdgesOf(i, [&](node i, node j, double value) {
		row[j] = value;
	});

	return row;
}

Vector Matrix::column(const index &j) const {
	Vector column(numberOfRows());
	graph.parallelForNodes([&](node i) {
		column[i] = graph.weight(i, j);
	});

	return column;
}

Matrix Matrix::operator+(const Matrix &other) const {
	return Matrix(*this) += other;
}

Matrix& Matrix::operator+=(const Matrix &other) {
	if (numberOfRows() != other.numberOfRows() || numberOfColumns() != other.numberOfColumns()) {
		throw std::runtime_error("dimensions of matrices do not match");
	}

	other.parallelForNonZeroElementsInRowOrder([&](node i, node j, double value) {
		if (i >= j) { // symmetric matrix
			graph.increaseWeight(i, j, value);
		}
	});

	return *this;
}

Matrix Matrix::operator-(const Matrix &other) const {
	return Matrix(*this) -= other;
}

Matrix& Matrix::operator-=(const Matrix &other) {
	if (numberOfRows() != other.numberOfRows() || numberOfColumns() != other.numberOfColumns()) {
		throw std::runtime_error("dimensions of matrices do not match");
	}

	other.parallelForNonZeroElementsInRowOrder([&](node i, node j, double value) {
		if (i >= j) { // symmetric matrix
			graph.increaseWeight(i, j, -value);
		}
	});

	return *this;
}

Matrix Matrix::operator*(const double &scalar) const {
	return Matrix(*this) *= scalar;
}

Matrix& Matrix::operator*=(const double &scalar) {
	graph.parallelForWeightedEdges([&](node i, node j, double value) {
		graph.setWeight(i, j, value * scalar);
	});

	return *this;
}

Vector Matrix::operator*(const Vector &vector) const {
	if (vector.isTransposed() && numberOfColumns() != 1) {
		throw std::runtime_error("vector is not transposed correctly");
	} else if (numberOfColumns() != vector.getDimension()) {
		throw std::runtime_error("dimensions of matrix and vector do not match");
	}

	Vector result(numberOfRows(), 0.0);

	parallelForNonZeroElementsInRowOrder([&](node i, node j, double value) {
		result[i] += value * vector[j];
	});

	return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
	if (numberOfColumns() != other.numberOfRows()) {
		throw std::runtime_error("dimensions of matrices do not match");
	}

	std::vector<Vector> columns(other.numberOfColumns());
#pragma omp parallel for
	for (count j = 0; j < other.numberOfColumns(); ++j) {
		columns[j] = (*this * other.column(j));
	}

	return Matrix(columns);
}

Matrix Matrix::operator/(const double &divisor) const {
	return Matrix(*this) /= divisor;
}

Matrix& Matrix::operator/=(const double &divisor) {
	return *this *= 1 / divisor;
}



} /* namespace NetworKit */
