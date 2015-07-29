/*
 * CSRMatrix.h
 *
 *  Created on: May 6, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef CSRMATRIX_H_
#define CSRMATRIX_H_

#include <vector>
#include "../Globals.h"
#include "Vector.h"
#include "../graph/Graph.h"
#include "../algebraic/SparseAccumulator.h"

namespace NetworKit {

class CSRMatrix {
private:
	std::vector<index> rowIdx;
	std::vector<index> columnIdx;
	std::vector<double> nonZeros;

	count nRows;
	count nCols;

public:
	/** Default constructor */
	CSRMatrix();

	CSRMatrix(const count nRows, const count nCols, const std::vector<std::pair<index, index>> &positions, const std::vector<double> &values);

	CSRMatrix(const count nRows, const count nCols, const std::vector<std::vector<index>> &columnIdx, const std::vector<std::vector<double>> &values);

	CSRMatrix (const CSRMatrix &other) = default;

	CSRMatrix (CSRMatrix &&other) = default;

	virtual ~CSRMatrix() = default;

	CSRMatrix& operator=(CSRMatrix &&other) = default;

	CSRMatrix& operator=(const CSRMatrix &other) = default;

	/**
	 * @return Number of rows.
	 */
	inline count numberOfRows() const {
		return nRows;
	}

	/**
	 * @return Number of columns.
	 */
	inline count numberOfColumns() const {
		return nCols;
	}

	/**
	 * @param i The row index.
	 * @return Number of non-zeros in row @a i.
	 */
	count nnzInRow(const index i) const;

	/**
	 * @return Number of non-zeros in this matrix.
	 */
	count nnz() const;

	/**
	 * @return Value at matrix position (i,j).
	 */
	double operator()(const index i, const index j) const;

	/**
	 * Set the matrix at position (@a i, @a j) to @a value.
	 */
	void setValue(const index i, const index j, const double value);

	/**
	 * @return Row @a i of this matrix as vector.
	 */
	Vector row(const index i) const;

	/**
	 * @return Column @a j of this matrix as vector.
	 */
	Vector column(const index j) const;

	/**
	 * @return The main diagonal of this matrix.
	 */
	Vector diagonal() const;

	/**
	 * Adds this matrix to @a other and returns the result.
	 * @return The sum of this matrix and @a other.
	 */
	CSRMatrix operator+(const CSRMatrix &other) const;

	/**
	 * Adds @a other to this matrix.
	 * @return Reference to this matrix.
	 */
	CSRMatrix& operator+=(const CSRMatrix &other);

	/**
	 * Subtracts @a other from this matrix and returns the result.
	 * @return The difference of this matrix and @a other.
	 *
	 */
	CSRMatrix operator-(const CSRMatrix &other) const;

	/**
	 * Subtracts @a other from this matrix.
	 * @return Reference to this matrix.
	 */
	CSRMatrix& operator-=(const CSRMatrix &other);

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar and returns the result.
	 * @return The result of multiplying this matrix with @a scalar.
	 */
	CSRMatrix operator*(const double &scalar) const;

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar.
	 * @return Reference to this matrix.
	 */
	CSRMatrix& operator*=(const double &scalar);

	/**
	 * Multiplies this matrix with @a vector and returns the result.
	 * @return The result of multiplying this matrix with @a vector.
	 */
	Vector operator*(const Vector &vector) const;

	/**
	 * Multiplies this matrix with @a other and returns the result in a new matrix.
	 * @return The result of multiplying this matrix with @a other.
	 */
	CSRMatrix operator*(const CSRMatrix &other) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor and returns the result in a new matrix.
	 * @return The result of dividing this matrix by @a divisor.
	 */
	CSRMatrix operator/(const double &divisor) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor.
	 * @return Reference to this matrix.
	 */
	CSRMatrix& operator/=(const double &divisor);

	template<typename L> static CSRMatrix binaryOperator(const CSRMatrix &A, const CSRMatrix &B, L binaryOp);

	static CSRMatrix mTmMultiply(const CSRMatrix &A, const CSRMatrix &B);

	static CSRMatrix mmTMultiply(const CSRMatrix &A, const CSRMatrix &B);

	static Vector mTvMultiply(const CSRMatrix &matrix, const Vector &vector);

	static CSRMatrix graphLaplacian(const Graph &graph);

	static CSRMatrix adjacencyMatrix(const Graph &graph);

	static Graph laplacianToGraph(const CSRMatrix &laplacian);

	CSRMatrix transpose() const;

	/**
	 * Iterate over all non-zero elements of row @a row in the matrix and call handler(index row, index column, double value)
	 */
	template<typename L> void forNonZeroElementsInRow(index i, L handle) const;

	template<typename L> void parallelForNonZeroElementsInRow(index i, L handle) const;

	/**
	 * Iterate over all non-zero elements of the matrix in row order and call handler (lambda closure).
	 */
	template<typename L> void forNonZeroElementsInRowOrder(L handle) const;

	/**
	 * Iterate in parallel over all rows and call handler (lambda closure) on non-zero elements of the matrix.
	 */
	template<typename L> void parallelForNonZeroElementsInRowOrder(L handle) const;

	/**
	 * Iterate in parallel over all rows and call handler (lambda closure) on non-zero elements of the matrix.
	 */
	template<typename L> void parallelForNonZeroElementsInRowOrder(L handle);
};

template<typename L> inline CSRMatrix NetworKit::CSRMatrix::binaryOperator(const CSRMatrix &A, const CSRMatrix &B, L binaryOp) {
	assert(A.nRows == B.nRows && A.nCols == B.nCols);
	std::vector<int64_t> columnPointer(A.nCols, -1);
	std::vector<double> Arow(A.nCols, 0.0);
	std::vector<double> Brow(A.nCols, 0.0);

	std::vector<std::pair<index,index>> positions;
	std::vector<double> values;

	for (index i = 0; i < A.nRows; ++i) {
		index listHead = 0;
		count nnz = 0;

		// search for nonZeros in our own matrix
		for (index k = A.rowIdx[i]; k < A.rowIdx[i+1]; ++k) {
			index j = A.columnIdx[k];
			Arow[j] = A.nonZeros[k];

			columnPointer[j] = listHead;
			listHead = j;
			nnz++;
		}

		// search for nonZeros in the other matrix
		for (index k = B.rowIdx[i]; k < B.rowIdx[i+1]; ++k) {
			index j = B.columnIdx[k];
			Brow[j] = B.nonZeros[k];

			if (columnPointer[j] == -1) { // our own matrix does not have a nonZero entry in column j
				columnPointer[j] = listHead;
				listHead = j;
				nnz++;
			}
		}

		// apply operator on the found nonZeros in A and B
		for (count k = 0; k < nnz; ++k) {
			double value = binaryOp(Arow[listHead], Brow[listHead]);
			if (value != 0.0) {
				positions.push_back(std::make_pair(i, listHead));
				values.push_back(value);
			}

			index temp = listHead;
			listHead = columnPointer[listHead];

			// reset for next row
			columnPointer[temp] = -1;
			Arow[temp] = 0.0;
			Brow[temp] = 0.0;
		}

		nnz = 0;
	}

	return CSRMatrix(A.nRows, A.nCols, positions, values);
}

} /* namespace NetworKit */

template<typename L>
inline void NetworKit::CSRMatrix::forNonZeroElementsInRow(index i, L handle) const {
	for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
		handle(columnIdx[k], nonZeros[k]);
	}
}

template<typename L>
inline void NetworKit::CSRMatrix::parallelForNonZeroElementsInRow(index i, L handle) const {
#pragma omp parallel for
	for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
		handle(columnIdx[k], nonZeros[k]);
	}
}

template<typename L>
inline void NetworKit::CSRMatrix::forNonZeroElementsInRowOrder(L handle) const {
	for (index i = 0; i < nRows; ++i) {
		for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
			handle(i, columnIdx[k], nonZeros[k]);
		}
	}
}

template<typename L>
inline void NetworKit::CSRMatrix::parallelForNonZeroElementsInRowOrder(L handle) const {
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
			handle(i, columnIdx[k], nonZeros[k]);
		}
	}
}

template<typename L>
inline void NetworKit::CSRMatrix::parallelForNonZeroElementsInRowOrder(L handle) {
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		for (index k = rowIdx[i]; k < rowIdx[i+1]; ++k) {
			handle(i, columnIdx[k], nonZeros[k]);
		}
	}
}

#endif /* TESTMATRIX_H_ */
