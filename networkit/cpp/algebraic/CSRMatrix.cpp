/*
 * CSRMatrix.cpp
 *
 *  Created on: May 6, 2015
 *      Author: Michael
 */

#include "CSRMatrix.h"
#include "SparseAccumulator.h"

#include <cassert>

namespace NetworKit {

CSRMatrix::CSRMatrix() : rowIdx(0), columnIdx(0), nonZeros(0), nRows(0), nCols(0) {
}

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const std::vector<std::pair<index, index>> &positions, const std::vector<double> &values) : nRows(nRows), nCols(nCols) {
	count nnz = values.size();
	rowIdx = std::vector<index>(nRows + 1, 0);
	columnIdx = std::vector<index>(nnz);
	nonZeros = std::vector<double>(nnz);

	for (index i = 0; i < nnz; ++i) {
		rowIdx[positions[i].first]++;
	}

	for (index i = 0, prefixSum = 0; i < nRows; ++i) {
		count nnzInRow = rowIdx[i];
		rowIdx[i] = prefixSum;
		prefixSum += nnzInRow;
	}
	rowIdx[nRows] = nnz;

	for (index i = 0; i < nnz; ++i) {
		index row = positions[i].first;
		index dest = rowIdx[row];

		columnIdx[dest] = positions[i].second;
		nonZeros[dest] = values[i];

		rowIdx[row]++;
	}

	for (index i = 0, firstIdxOfRow = 0; i <= nRows; ++i) {
		index newRow = rowIdx[i];
		rowIdx[i] = firstIdxOfRow;
		firstIdxOfRow = newRow;
	}
}

CSRMatrix::CSRMatrix(const count nRows, const count nCols, const std::vector<std::vector<index>> &columnIdx, const std::vector<std::vector<double>> &values) : nRows(nRows), nCols(nCols) {
	 rowIdx = std::vector<index>(nRows + 1, 0);
	 count nnz = columnIdx[0].size();
	 for (index i = 1; i < columnIdx.size(); ++i) {
		 rowIdx[i] = rowIdx[i-1] + columnIdx[i-1].size();
		 nnz += columnIdx[i].size();
	 }
	 rowIdx[nRows] = nnz;

	 this->columnIdx = std::vector<index>(nnz);
	 this->nonZeros = std::vector<double>(nnz);

#pragma omp parallel for
	 for (index i = 0; i < nRows; ++i) {
		 for (index k = 0; k < columnIdx[i].size(); ++k) {
			 this->columnIdx[rowIdx[i] + k] = columnIdx[i][k];
			 nonZeros[rowIdx[i] + k] = values[i][k];
		 }
	 }
}

count CSRMatrix::nnzInRow(const index i) const {
	assert(i >= 0 && i < nRows);
	return rowIdx[i+1] - rowIdx[i];
}

count CSRMatrix::nnz() const {
	return nonZeros.size();
}

double CSRMatrix::operator()(const index i, const index j) const {
	assert(i >= 0 && i < nRows);
	assert(j >= 0 && j < nCols);

	double value = 0.0;
	for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
		if (columnIdx[k] == j) {
			value = nonZeros[k];
			break;
		}
	}

	return value;
}

void CSRMatrix::setValue(const index i, const index j, const double value) {
//	assert(i >= 0 && i < nRows);
//	assert(j >= 0 && j < nCols);
//
//
//
//	graph.setWeight(i, j, value);
}

Vector CSRMatrix::row(const index i) const {
	assert(i >= 0 && i < nRows);

	Vector row(numberOfColumns(), 0.0, true);
	parallelForNonZeroElementsInRow(i, [&](index j, double value) {
		row[j] = value;
	});

	return row;
}

Vector CSRMatrix::column(const index j) const {
	assert(j >= 0 && j < nCols);

	Vector column(numberOfRows());
#pragma omp parallel for
	for (node i = 0; i < numberOfRows(); ++i) {
		column[i] = (*this)(i,j);
	}

	return column;
}

Vector CSRMatrix::diagonal() const {
	Vector diag(std::min(nRows, nCols), 0.0);
	for (index i = 0; i < diag.getDimension(); ++i) {
		diag[i] = (*this)(i,i);
	}

	return diag;
}

CSRMatrix CSRMatrix::operator+(const CSRMatrix &other) const {
	assert(nRows == other.nRows && nCols == other.nCols);
	return CSRMatrix::binaryOperator(*this, other, [&](double val1, double val2) {return val1 + val2;});
}

CSRMatrix& CSRMatrix::operator+=(const CSRMatrix &other) {
	assert(nRows == other.nRows && nCols == other.nCols);
	*this = CSRMatrix::binaryOperator(*this, other, [&](double val1, double val2) {return val1 + val2;});
	return *this;
}

CSRMatrix CSRMatrix::operator-(const CSRMatrix &other) const {
	assert(nRows == other.nRows && nCols == other.nCols);
	return CSRMatrix::binaryOperator(*this, other, [&](double val1, double val2) {return val1 - val2;});
}

CSRMatrix& CSRMatrix::operator-=(const CSRMatrix &other) {
	assert(nRows == other.nRows && nCols == other.nCols);
	*this = CSRMatrix::binaryOperator(*this, other, [&](double val1, double val2) {return val1 + val2;});
	return *this;
}

CSRMatrix CSRMatrix::operator*(const double &scalar) const {
	return CSRMatrix(*this) *= scalar;
}

CSRMatrix& CSRMatrix::operator*=(const double &scalar) {
#pragma omp parallel for
	for (index k = 0; k < nonZeros.size(); ++k) {
		nonZeros[k] *= scalar;
	}

	return *this;
}

Vector CSRMatrix::operator*(const Vector &vector) const {
	assert(!vector.isTransposed());
	assert(nCols == vector.getDimension());

	Vector result(nRows, 0.0);
	parallelForNonZeroElementsInRowOrder([&](node i, node j, double value) {
		result[i] += value * vector[j];
	});

	return result;
}

CSRMatrix CSRMatrix::operator*(const CSRMatrix &other) const {
	assert(nCols == other.nRows);

	std::vector<std::pair<index,index>> positions;
	std::vector<double> values;

	SparseAccumulator spa(numberOfRows());
	for (index i = 0; i < numberOfRows(); ++i) {
		forNonZeroElementsInRow(i, [&](index k, double val1) {
			other.forNonZeroElementsInRow(k, [&](index j, double val2) {
				spa.scatter(val1 * val2, j);
			});
		});

		spa.gather([&](index i, index j, double value){
			positions.push_back(std::make_pair(i,j));
			values.push_back(value);
		});

		spa.increaseRow();
	}

	return CSRMatrix(nRows, other.nCols, positions, values);
}

CSRMatrix CSRMatrix::operator/(const double &divisor) const {
	return CSRMatrix(*this) /= divisor;
}

CSRMatrix& CSRMatrix::operator/=(const double &divisor) {
	return *this *= 1.0 / divisor;
}



CSRMatrix CSRMatrix::mTmMultiply(const CSRMatrix &A, const CSRMatrix &B) {
	assert(A.nRows == B.nRows);

	std::vector<std::vector<index>> columnIdx(A.numberOfColumns());
	std::vector<std::vector<double>> values(A.numberOfColumns());

	for (index k = 0; k < A.numberOfRows(); ++k) {
		A.forNonZeroElementsInRow(k, [&](index i, double vA) {
			B.forNonZeroElementsInRow(k, [&](index j, double vB) {
				bool found = false;
				for (index l = 0; l < columnIdx[i].size(); ++l) {
					if (columnIdx[i][l] == j) {
						values[i][l] += vA * vB;
						found = true;
						break;
					}
				}

				if (!found) {
					columnIdx[i].push_back(j);
					values[i].push_back(vA * vB);
				}
			});
		});
	}

	return CSRMatrix(A.nCols, B.nCols, columnIdx, values);
}

CSRMatrix CSRMatrix::mmTMultiply(const CSRMatrix &A, const CSRMatrix &B) {
	assert(A.nCols == B.nCols);

	std::vector<std::vector<index>> columnIdx(A.numberOfRows());
	std::vector<std::vector<double>> values(A.numberOfRows());

	for (index i = 0; i < A.numberOfRows(); ++i) {
		A.forNonZeroElementsInRow(i, [&](index k, double vA) {
			for (index j = 0; j < B.numberOfRows(); ++j) {
				double vB = B(j,k);
				if (vB != 0.0) {
					bool found = false;
					for (index l = 0; l < columnIdx[i].size(); ++l) {
						if (columnIdx[i][l] == j) {
							values[i][l] += vA * vB;
							found = true;
							break;
						}
					}

					if (!found) {
						columnIdx[i].push_back(j);
						values[i].push_back(vA * vB);
					}
				}
			}
		});
	}

	return CSRMatrix(A.nRows, B.nRows, columnIdx, values);
}

Vector CSRMatrix::mTvMultiply(const CSRMatrix &matrix, const Vector &vector) {
	assert(matrix.nRows == vector.getDimension());

	Vector result(matrix.numberOfColumns(), 0.0);
	for (index k = 0; k < matrix.numberOfRows(); ++k) {
		matrix.forNonZeroElementsInRow(k, [&](index j, double value) {
			result[j] += value * vector[k];
		});
	}

	return result;
}

CSRMatrix CSRMatrix::graphLaplacian(const Graph &graph) {
	std::vector<std::pair<index, index>> positions;
	std::vector<double> values;

	graph.forNodes([&](const index i){
		double weightedDegree = 0.0;

		double selfLoopWeight = 0.0;
		graph.forNeighborsOf(i, [&](const index j, double weight) { // - adjacency matrix
			if (j == i) {
				selfLoopWeight = weight;
			} else {
				positions.push_back(std::make_pair(i,j));
				values.push_back(-weight);
			}

			weightedDegree += weight;
		});

		positions.push_back(std::make_pair(i,i));
		values.push_back(weightedDegree - selfLoopWeight); // degree matrix
	});

	return CSRMatrix(graph.upperNodeIdBound(), graph.upperNodeIdBound(), positions, values);
}

CSRMatrix CSRMatrix::adjacencyMatrix(const Graph &graph) {
	int nonZeros = graph.isDirected()? graph.numberOfEdges() : graph.numberOfEdges() * 2;

	std::vector<std::pair<index, index>> positions(nonZeros);
	std::vector<double> values(nonZeros);

	int index = 0;
	graph.forEdges([&](node i, node j, double val) {
		positions[index] = std::make_pair(i,j);
		values[index] = val;
		index++;
		if (!graph.isDirected() && i != j) {
			positions[index] = std::make_pair(j,i);
			values[index] = val;
			index++;
		}
	});

	return CSRMatrix(graph.numberOfNodes(), graph.numberOfNodes(), positions, values);
}

Graph CSRMatrix::laplacianToGraph(const CSRMatrix &laplacian) {
	Graph G(std::max(laplacian.numberOfRows(), laplacian.numberOfColumns()), true, true);
	laplacian.forNonZeroElementsInRowOrder([&](node u, node v, edgeweight weight) {
		G.addEdge(u, v, -weight);
	});

	return G;
}

CSRMatrix CSRMatrix::transpose() const {
//	Matrix transposedMatrix(numberOfColumns(), numberOfRows());
//	parallelForNonZeroElementsInRowOrder([&](index i, index j, edgeweight weight){
//		transposedMatrix.graph.addEdge(i,j,weight);
//	});
//
//	return transposedMatrix;
}



} /* namespace NetworKit */
