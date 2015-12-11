/*
 * DenseMatrix.h
 *
 *  Created on: Nov 25, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_DENSEMATRIX_H_
#define NETWORKIT_CPP_ALGEBRAIC_DENSEMATRIX_H_

#include "../Globals.h"
#include "Vector.h"

#include <cassert>
#include <vector>

namespace NetworKit {

/**
 * Represents a dense matrix. Use this matrix to run LU decompositions and LU solves.
 * Note that most matrices are rather sparse s.t. CSRMatrix might be a better representation.
 */
class DenseMatrix {
private:
	count nRows;
	count nCols;
	std::vector<double> entries;

public:
	/** Default constructor */
	DenseMatrix();

	/**
	 * Constructs an instance of DenseMatrix given the number of rows (@a nRows) and the number of columns (@a nCols) and its
	 * values (@a entries).
	 * @param nRows Number of rows.
	 * @param nCols Number of columns.
	 * @param entries Entries of the matrix.
	 * @note The size of the @a entries vector should be equal to @a nRows * @a nCols.
	 */
	DenseMatrix(const count nRows, const count nCols, const std::vector<double> &entries);

	/** Default destructor */
	virtual ~DenseMatrix() = default;

	/** Default copy constructor */
	DenseMatrix (const DenseMatrix &other) = default;

	/** Default move constructor */
	DenseMatrix (DenseMatrix &&other) = default;

	/** Default copy assignment operator */
	DenseMatrix& operator=(DenseMatrix &&other) = default;

	/** Default move assignment operator */
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

	/**
	 * Decomposes the given @a matrix into lower L and upper U matrix (in-place).
	 * @param matrix The matrix to decompose into LU.
	 */
	static void LUDecomposition(DenseMatrix &matrix);

	/**
	 * Computes the solution vector x to the system @a LU * x = @a b where @a LU is a matrix decomposed into L and U.
	 * @param LU Matrix decomposed into lower L and upper U matrix.
	 * @param b Right-hand side.
	 * @return Solution vector x to the linear equation system LU * x = b.
	 */
	static Vector LUSolve(const DenseMatrix &LU, const Vector &b);

	/**
	 * Computes @a A @a binaryOp @a B on the elements of matrix @a A and matrix @a B.
	 * @param A
	 * @param B
	 * @param binaryOp Function handling (double, double) -> double
	 * @return @a A @a binaryOp @a B.
	 * @note @a A and @a B must have the same dimensions.
	 */
	template<typename L> static DenseMatrix binaryOperator(const DenseMatrix &A, const DenseMatrix &B, L binaryOp);

	/**
	 * Iterate over all non-zero elements of row @a row in the matrix and call handler(index column, double value)
	 */
	template<typename L> void forElementsInRow(index i, L handle) const;

	/**
	 * Iterate in parallel over all non-zero elements of row @a row in the matrix and call handler(index column, double value)
	 */
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
