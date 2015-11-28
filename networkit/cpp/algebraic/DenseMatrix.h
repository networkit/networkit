/*
 * DenseMatrix.h
 *
 *  Created on: Nov 25, 2015
 *      Author: Michael
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_DENSEMATRIX_H_
#define NETWORKIT_CPP_ALGEBRAIC_DENSEMATRIX_H_

#include "../Globals.h"
#include "Vector.h"

#include <cassert>
#include <vector>

namespace NetworKit {

class DenseMatrix {
private:
	count nRows;
	count nCols;
	std::vector<double> entries;

public:
	/** Default constructor */
	DenseMatrix();

	DenseMatrix(const count nRows, const count nCols, const std::vector<double> &entries);

	virtual ~DenseMatrix() = default;

	DenseMatrix (const DenseMatrix &other) = default;

	DenseMatrix (DenseMatrix &&other) = default;

	DenseMatrix& operator=(DenseMatrix &&other) = default;

	DenseMatrix& operator=(const DenseMatrix &other) = default;

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
	DenseMatrix operator+(const DenseMatrix &other) const;

	/**
	 * Adds @a other to this matrix.
	 * @return Reference to this matrix.
	 */
	DenseMatrix& operator+=(const DenseMatrix &other);

	/**
	 * Subtracts @a other from this matrix and returns the result.
	 * @return The difference of this matrix and @a other.
	 *
	 */
	DenseMatrix operator-(const DenseMatrix &other) const;

	/**
	 * Subtracts @a other from this matrix.
	 * @return Reference to this matrix.
	 */
	DenseMatrix& operator-=(const DenseMatrix &other);

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar and returns the result.
	 * @return The result of multiplying this matrix with @a scalar.
	 */
	DenseMatrix operator*(const double &scalar) const;

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar.
	 * @return Reference to this matrix.
	 */
	DenseMatrix& operator*=(const double &scalar);

	/**
	 * Multiplies this matrix with @a vector and returns the result.
	 * @return The result of multiplying this matrix with @a vector.
	 */
	Vector operator*(const Vector &vector) const;

	/**
	 * Multiplies this matrix with @a other and returns the result in a new matrix.
	 * @return The result of multiplying this matrix with @a other.
	 */
	DenseMatrix operator*(const DenseMatrix &other) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor and returns the result in a new matrix.
	 * @return The result of dividing this matrix by @a divisor.
	 */
	DenseMatrix operator/(const double &divisor) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor.
	 * @return Reference to this matrix.
	 */
	DenseMatrix& operator/=(const double &divisor);

	static void LUDecomposition(DenseMatrix &matrix);

	static Vector LUSolve(const DenseMatrix &LU, const Vector &b);

	template<typename L> static DenseMatrix binaryOperator(const DenseMatrix &A, const DenseMatrix &B, L binaryOp);

	/**
	 * Iterate over all non-zero elements of row @a row in the matrix and call handler(index row, index column, double value)
	 */
	template<typename L> void forElementsInRow(index i, L handle) const;

	template<typename L> void parallelForElementsInRow(index i, L handle) const;

	/**
	 * Iterate over all non-zero elements of the matrix in row order and call handler (lambda closure).
	 */
	template<typename L> void forElementsInRowOrder(L handle) const;

	/**
	 * Iterate in parallel over all rows and call handler (lambda closure) on non-zero elements of the matrix.
	 */
	template<typename L> void parallelForElementsInRowOrder(L handle) const;

	/**
	 * Iterate in parallel over all rows and call handler (lambda closure) on non-zero elements of the matrix.
	 */
	template<typename L> void parallelForElementsInRowOrder(L handle);
};

template<typename L> inline DenseMatrix NetworKit::DenseMatrix::binaryOperator(const DenseMatrix &A, const DenseMatrix &B, L binaryOp) {
	assert(A.nRows == B.nRows && A.nCols == B.nCols);

	std::vector<double> resultEntries(A.numberOfRows() * A.numberOfColumns());

#pragma omp parallel for
	for (index i = 0; i < A.numberOfRows(); ++i) {
		index offset = i * A.numberOfColumns();
		for (index j = offset; j < offset + A.numberOfColumns(); ++j) {
			resultEntries[j] = binaryOp(A.entries[j], B.entries[j]);
		}
	}

	return DenseMatrix(A.numberOfRows(), A.numberOfColumns(), resultEntries);
}

template<typename L>
inline void NetworKit::DenseMatrix::forElementsInRow(index i, L handle) const {
	index offset = i * numberOfColumns();
	for (index k = offset, j = 0; k < offset + numberOfColumns(); ++k, ++j) {
		handle(j, entries[k]);
	}
}

template<typename L>
inline void NetworKit::DenseMatrix::parallelForElementsInRow(index i, L handle) const {
	index offset = i * numberOfColumns();
#pragma omp parallel for
	for (index j = 0; j < numberOfColumns(); ++j) {
		handle(j, entries[offset + j]);
	}
}

template<typename L>
inline void NetworKit::DenseMatrix::forElementsInRowOrder(L handle) const {
	for (index i = 0; i < nRows; ++i) {
		index offset = i * numberOfColumns();
		for (index k = offset, j = 0; k < offset + numberOfColumns(); ++k, ++j) {
			handle(i, j, entries[k]);
		}
	}
}

template<typename L>
inline void NetworKit::DenseMatrix::parallelForElementsInRowOrder(L handle) const {
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		index offset = i * numberOfColumns();
		for (index k = offset, j = 0; k < offset + numberOfColumns(); ++k, ++j) {
			handle(i, j, entries[k]);
		}
	}
}

template<typename L>
inline void NetworKit::DenseMatrix::parallelForElementsInRowOrder(L handle) {
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		index offset = i * numberOfColumns();
		for (index k = offset, j = 0; k < offset + numberOfColumns(); ++k, ++j) {
			handle(i, j, entries[k]);
		}
	}
}

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_ALGEBRAIC_DENSEMATRIX_H_ */
