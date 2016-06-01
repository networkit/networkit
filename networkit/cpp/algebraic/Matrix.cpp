/*
 * Matrix.cpp
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "Matrix.h"

namespace NetworKit {

Matrix::Matrix() : graph(0, true, true), nRows(0), nCols(0), zero(0.0) {
}

Matrix::Matrix(const count dimension, const double zero) : graph(dimension, true, true), nRows(dimension), nCols(dimension), zero(zero) {
}

Matrix::Matrix(const count nRows, const count nCols, const double zero) : graph(std::max(nRows, nCols), true, true), nRows(nRows), nCols(nCols), zero(zero) {
}

Matrix::Matrix(const count dimension, const std::vector<Triplet>& triplets, const double zero) : Matrix(dimension, dimension, triplets, zero) {
}

Matrix::Matrix(const count nRows, const count nCols, const std::vector<Triplet>& triplets, const double zero) : graph(std::max(nRows, nCols), true, true), nRows(nRows), nCols(nCols), zero(zero) {
	for (size_t k = 0; k < triplets.size(); ++k) {
		assert(triplets[k].row < nRows && triplets[k].column < nCols);
		graph.addEdge(triplets[k].row, triplets[k].column, triplets[k].value);
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

	double weight = graph.weight(i,j);
	return weight == nullWeight? zero : weight;
}

void Matrix::setValue(const index i, const index j, const double value) {
	assert(i >= 0 && i < nRows);
	assert(j >= 0 && j < nCols);

	graph.setWeight(i, j, value);
}

Vector Matrix::row(const index i) const {
	assert(i >= 0 && i < nRows);

	Vector row(numberOfColumns(), zero, true);
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
		double val = graph.weight(i,j);
		column[i] = (val == nullWeight? zero : val);
	}

	return column;
}

Vector Matrix::diagonal() const {
	Vector diag(std::min(nRows, nCols), zero);
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
	Vector result(numberOfRows(), zero);

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
				if (wB != A.zero) {
					C.graph.increaseWeight(i, j, wA * wB);
				}
			}
		});
	}

	return C;
}

Vector Matrix::mTvMultiply(const Matrix &matrix, const Vector &vector) {
	assert(matrix.nRows == vector.getDimension());

	Vector result(matrix.numberOfColumns(), matrix.getZero());
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

Matrix Matrix::adjacencyMatrix(const Graph& graph) {
	Matrix A(graph.upperNodeIdBound());
	graph.forEdges([&](node u, node v, edgeweight w) {
		A.setValue(u, v, w);
		if (!graph.isDirected()) { // add symmetric value at (v, u)
			A.setValue(v, u, w);
		}
	});

	return A;
}

Matrix Matrix::diagonalMatrix(const Vector& diagonalElements) {
	Matrix D(diagonalElements.getDimension());
	for (index i = 0; i < diagonalElements.getDimension(); ++i) {
		D.setValue(i, i, diagonalElements[i]);
	}

	return D;
}

Matrix Matrix::incidenceMatrix(const Graph& graph) {
	if (!graph.hasEdgeIds()) throw std::runtime_error("Graph has no edge Ids. Index edges first by calling graph.indexEdges()");
	Matrix I(graph.upperNodeIdBound(), graph.upperEdgeIdBound());
	if (graph.isDirected()) {
		graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId) {
			if (u != v) {
				edgeweight w = sqrt(weight);
				I.setValue(u, edgeId, w);
				I.setValue(v, edgeId, -w);
			}
		});
	} else {
		graph.forEdges([&](node u, node v, edgeweight weight, edgeid edgeId){
			if (u != v) {
				edgeweight w = sqrt(weight);
				if (u < v) { // orientation: small node number -> great node number
					I.setValue(u, edgeId, w);
					I.setValue(v, edgeId, -w);
				} else {
					I.setValue(u, edgeId, -w);
					I.setValue(v, edgeId, w);
				}
			}
		});
	}

	return I;
}

Matrix Matrix::laplacianMatrix(const Graph& graph) {
	Matrix L(graph.upperNodeIdBound());
	graph.forNodes([&](const index i){
		double weightedDegree = 0.0;

		graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
			L.setValue(i, j, -weight);
			if (i != j) { // exclude weight of diagonal since it would be subtracted later
				weightedDegree += weight;
			}
		});

		L.setValue(i, i, weightedDegree); // degree matrix
	});

	return L;
}

Matrix Matrix::normalizedLaplacianMatrix(const Graph& graph) {
	Matrix nL(graph.upperNodeIdBound());

	std::vector<double> weightedDegrees(graph.upperNodeIdBound(), 0.0);
	graph.parallelForNodes([&](const node u) {
		weightedDegrees[u] = graph.weightedDegree(u);
	});

	graph.forNodes([&](const node i){
		graph.forNeighborsOf(i, [&](const node j, double weight){
			if (i != j) {
				nL.setValue(i, j, -weight/sqrt(weightedDegrees[i] * weightedDegrees[j]));
			}
		});

		if (weightedDegrees[i] != 0.0) {
			if (graph.isWeighted()) {
				nL.setValue(i, i, 1-(graph.weight(i, i)) / weightedDegrees[i]);
			} else {
				nL.setValue(i, i, 1);
			}
		}
	});

	return nL;
}





} /* namespace NetworKit */
