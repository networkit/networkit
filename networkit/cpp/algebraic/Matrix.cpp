/*
 * Matrix.cpp
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "Matrix.h"

namespace NetworKit {

Matrix::Matrix() : graph(0, true, true), nRows(0), nCols(0) {
}

Matrix::Matrix(const count &dimension) : graph(dimension, true, true), nRows(dimension), nCols(dimension) {
}

Matrix::Matrix(const count &nRows, const count &nCols) : graph(std::max(nRows, nCols), true, true), nRows(nRows), nCols(nCols) {
}

Matrix::Matrix(const count &dimension, const std::vector<std::pair<index, index>> &positions,
					const std::vector<double> &values) : graph(dimension, true, true), nRows(dimension), nCols(dimension) {
	assert(positions.size() == values.size());

	for (size_t k = 0; k < positions.size(); ++k) {
		assert(positions[k].first >= 0 && positions[k].second >= 0 && positions[k].first < dimension && positions[k].second < dimension);

		std::pair<node, node> pos = positions[k];
		graph.addEdge(pos.first, pos.second, values[k]);
	}
}

Matrix::Matrix(const count &nRows, const count &nCols, const std::vector<std::pair<index, index>> &positions,
					const std::vector<double> &values) : graph(std::max(nRows, nCols), true, true), nRows(nRows), nCols(nCols) {
	assert(positions.size() == values.size());

	for (size_t k = 0; k < positions.size(); ++k) {
		assert(positions[k].first >= 0 && positions[k].second >= 0 && positions[k].first < nRows && positions[k].second < nCols);

		std::pair<node, node> pos = positions[k];
		graph.addEdge(pos.first, pos.second, values[k]);
	}
}

Matrix::Matrix(const std::vector<Vector> &rows) {
	assert(rows.size() > 0);

	nRows = rows.size();
	nCols = rows[0].getDimension();
	graph = Graph(std::max(nRows, nCols), true, true);

#pragma omp parallel for
	for (size_t i = 0; i < nRows; ++i) {
		if (rows[i].getDimension() != nCols) {
			throw std::runtime_error("Matrix(const std::vector<Vector> &rows): Column dimensions of one or more rows do not match");
		}
	}

	for (index i = 0; i < nRows; ++i) {
		for (index j = 0; j < nCols; ++j) {
			double value = rows[i][j];
			if (value != 0.0) { // do not store 0 values
				graph.addEdge(i, j, value);
			}
		}
	}
}

Matrix::Matrix(const Matrix &other) : graph(other.graph), nRows(other.nRows), nCols(other.nCols) {
}

count Matrix::numberOfRows() const {
	return nRows;
}

count Matrix::numberOfColumns() const {
	return nCols;
}

double Matrix::operator()(const index &i, const index &j) const {
	if (i < 0 || i >= numberOfRows()) {
		throw std::out_of_range("Matrix(i,j): Row index out of range");
	} else if (j< 0 || j >= numberOfColumns()) {
		throw std::out_of_range("Matrix(i,j): Column index out of range");
	}

	return graph.weight(i,j);
}

void Matrix::setValue(const index &i, const index &j, const double &value) {
	if (i < 0 || i >= numberOfRows()) {
		throw std::out_of_range("Matrix::setValue(const index &i, const index &j, const double &value): "
																						"Row index out of range");
	} else if (j< 0 || j >= numberOfColumns()) {
		throw std::out_of_range("Matrix::setValue(const index &i, const index &j, const double &value): "
																						"Column index out of range");
	}

	graph.setWeight(i, j, value);
}

Vector Matrix::row(const index &i) const {
	if (i < 0 || i >= numberOfRows()) {
		throw std::out_of_range("Matrix::row(const index &i): Row index out of range");
	}

	Vector row(numberOfColumns(), 0.0, true);
	graph.forEdgesOf(i, [&](node i, node j, double value) {
		row[j] = value;
	});

	return row;
}

Vector Matrix::column(const index &j) const {
	if (j < 0 || j >= numberOfColumns()) {
		throw std::out_of_range("Matrix::column(const index &j): Column index out of range");
	}

	Vector column(numberOfRows());
#pragma omp parallel for
	for (node i = 0; i < numberOfRows(); ++i) {
		column[i] = graph.weight(i,j);
	}

	return column;
}

Matrix Matrix::operator+(const Matrix &other) const {
	return Matrix(*this) += other;
}

Matrix& Matrix::operator+=(const Matrix &other) {
	if (numberOfRows() != other.numberOfRows() || numberOfColumns() != other.numberOfColumns()) {
		throw std::runtime_error("Matrix::operator+=(const Matrix &other): Dimensions of matrices do not match");
	}

	other.parallelForNonZeroElementsInRowOrder([&](node i, node j, double value) {
		graph.increaseWeight(i, j, value);
	});

	return *this;
}

Matrix Matrix::operator-(const Matrix &other) const {
	return Matrix(*this) -= other;
}

Matrix& Matrix::operator-=(const Matrix &other) {
	if (numberOfRows() != other.numberOfRows() || numberOfColumns() != other.numberOfColumns()) {
		throw std::runtime_error("Matrix::operator-=(const Matrix &other): Dimensions of matrices do not match");
	}

	other.parallelForNonZeroElementsInRowOrder([&](node i, node j, double value) {
		graph.increaseWeight(i, j, -value);
	});

	return *this;
}

Matrix Matrix::operator*(const double &scalar) const {
	return Matrix(*this) *= scalar;
}

Matrix& Matrix::operator*=(const double &scalar) {
	graph.parallelForEdges([&](node i, node j, double value) {
		graph.setWeight(i, j, value * scalar);
	});

	return *this;
}

Vector Matrix::operator*(const Vector &vector) const {
	if (vector.isTransposed() && numberOfColumns() != 1) {
		throw std::runtime_error("operator*(const Vector &vector): Vector is not transposed correctly");
	} else if (numberOfColumns() != vector.getDimension()) {
		throw std::runtime_error("operator*(const Vector &vector): Dimensions of matrix and vector do not match");
	}

	Vector result(numberOfRows(), 0.0);

	parallelForNonZeroElementsInRowOrder([&](node i, node j, double value) {
		result[i] += value * vector[j];
	});

	return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
	if (numberOfColumns() != other.numberOfRows()) {
		throw std::runtime_error("Matrix::operator*(const Matrix &other): Dimensions of matrices do not match");
	}

	Matrix result(numberOfRows(), other.numberOfColumns());
	SparseAccumulator spa(numberOfRows());
	for (index r = 0; r < numberOfRows(); ++r) {
		graph.forNeighborsOf(r, [&](node v, double w1){
			other.graph.forNeighborsOf(v, [&](node u, double w2){
				double value = w1 * w2;
				spa.scatter(value, u);
			});
		});

		spa.gather([&](node row, node column, double value){
			result.graph.addEdge(row, column, value);
		});

		spa.increaseRow();
	}

	return result;
}

Matrix Matrix::operator/(const double &divisor) const {
	return Matrix(*this) /= divisor;
}

Matrix& Matrix::operator/=(const double &divisor) {
	return *this *= 1 / divisor;
}



} /* namespace NetworKit */
