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

Matrix::Matrix(const count dimension) : graph(dimension, true, true), nRows(dimension), nCols(dimension) {
}

Matrix::Matrix(const count nRows, const count nCols) : graph(std::max(nRows, nCols), true, true), nRows(nRows), nCols(nCols) {
}

Matrix::Matrix(const count dimension, const std::vector<std::pair<index, index>> &positions,
					const std::vector<double> &values) : graph(dimension, true, true), nRows(dimension), nCols(dimension) {
	assert(positions.size() == values.size());

	for (size_t k = 0; k < positions.size(); ++k) {
		assert(positions[k].first >= 0 && positions[k].second >= 0 && positions[k].first < dimension && positions[k].second < dimension);

		std::pair<node, node> pos = positions[k];
		graph.addEdge(pos.first, pos.second, values[k]);
	}
}

Matrix::Matrix(const count nRows, const count nCols, const std::vector<std::pair<index, index>> &positions,
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

count Matrix::nnzInRow(const index i) const {
	assert(i >= 0 && i < nRows);
	return graph.degree(i);
}

count Matrix::nnz() const {
	count nnz = 0;
	for (index i = 0; i < nRows; ++i) {
		nnz += nnzInRow(i);
	}

	return nnz;
}

double Matrix::operator()(const index i, const index j) const {
	assert(i >= 0 && i < nRows);
	assert(j >= 0 && j < nCols);

	return graph.weight(i,j);
}

void Matrix::setValue(const index i, const index j, const double value) {
	assert(i >= 0 && i < nRows);
	assert(j >= 0 && j < nCols);

	graph.setWeight(i, j, value);
}

Vector Matrix::row(const index i) const {
	assert(i >= 0 && i < nRows);

	Vector row(numberOfColumns(), 0.0, true);
	graph.forEdgesOf(i, [&](node i, node j, double value) {
		row[j] = value;
	});

	return row;
}

Vector Matrix::column(const index j) const {
	assert(j >= 0 && j < nCols);

	Vector column(numberOfRows());
#pragma omp parallel for
	for (node i = 0; i < numberOfRows(); ++i) {
		column[i] = graph.weight(i,j);
	}

	return column;
}

Vector Matrix::diagonal() const {
	Vector diag(std::min(nRows, nCols), 0);
	for (index i = 0; i < diag.getDimension(); ++i) {
		diag[i] = (*this)(i,i);
	}

	return diag;
}

Matrix Matrix::operator+(const Matrix &other) const {
	return Matrix(*this) += other;
}

Matrix& Matrix::operator+=(const Matrix &other) {
	assert(nRows == other.nRows && nCols == other.nCols);

	other.forNonZeroElementsInRowOrder([&](node i, node j, double value) {
		graph.increaseWeight(i, j, value);
	});

	return *this;
}

Matrix Matrix::operator-(const Matrix &other) const {
	return Matrix(*this) -= other;
}

Matrix& Matrix::operator-=(const Matrix &other) {
	assert(nRows == other.nRows && nCols == other.nCols);

	other.forNonZeroElementsInRowOrder([&](node i, node j, double value) {
		graph.increaseWeight(i, j, -value);
	});

	return *this;
}

Matrix Matrix::operator*(const double scalar) const {
	return Matrix(*this) *= scalar;
}

Matrix& Matrix::operator*=(const double scalar) {
	graph.parallelForEdges([&](node i, node j, double value) {
		graph.setWeight(i, j, value * scalar);
	});

	return *this;
}

Vector Matrix::operator*(const Vector &vector) const {
	assert(!vector.isTransposed());
	assert(nCols == vector.getDimension());
	Vector result(numberOfRows(), 0.0);

	parallelForNonZeroElementsInRowOrder([&](node i, node j, double value) {
		result[i] += value * vector[j];
	});

	return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
	assert(nCols == other.nRows);

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

Matrix Matrix::operator/(const double divisor) const {
	return Matrix(*this) /= divisor;
}

Matrix& Matrix::operator/=(const double divisor) {
	return *this *= 1 / divisor;
}

Matrix Matrix::mTmMultiply(const Matrix &A, const Matrix &B) {
	assert(A.nRows == B.nRows);

	Matrix C(A.numberOfColumns(), B.numberOfColumns());
	for (index k = 0; k < A.numberOfRows(); ++k) {
		A.graph.forNeighborsOf(k, [&](index i, edgeweight wA) {
			B.graph.forNeighborsOf(k, [&](index j, edgeweight wB) {
				C.graph.increaseWeight(i, j, wA * wB);
			});
		});
	}

	return C;
}

Matrix Matrix::mmTMultiply(const Matrix &A, const Matrix &B) {
	assert(A.nCols == B.nCols);

	Matrix C(A.numberOfRows(), B.numberOfRows());
	for (index i = 0; i < A.numberOfRows(); ++i) {
		A.graph.forNeighborsOf(i, [&](index k, edgeweight wA){
			for (index j = 0; j < B.numberOfRows(); ++j) {
				edgeweight wB = B(j,k);
				if (wB != 0.0) {
					C.graph.increaseWeight(i, j, wA * wB);
				}
			}
		});
	}

	return C;
}

Vector Matrix::mTvMultiply(const Matrix &matrix, const Vector &vector) {
	assert(matrix.nRows == vector.getDimension());

	Vector result(matrix.numberOfColumns(), 0.0);
	for (index k = 0; k < matrix.numberOfRows(); ++k) {
		matrix.graph.forNeighborsOf(k, [&](index j, edgeweight w){
			result[j] += w * vector[k];
		});
	}

	return result;
}

Matrix Matrix::transpose() const {
	Matrix transposedMatrix(numberOfColumns(), numberOfRows());
	parallelForNonZeroElementsInRowOrder([&](index i, index j, edgeweight weight){
		transposedMatrix.graph.addEdge(i,j,weight);
	});

	return transposedMatrix;
}



} /* namespace NetworKit */
