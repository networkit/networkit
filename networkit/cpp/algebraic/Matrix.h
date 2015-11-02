/*
 * Matrix.h
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "../graph/Graph.h"
#include "Vector.h"
#include "SparseAccumulator.h"

namespace NetworKit {

/**
 * @ingroup algebraic
 * The matrix class represents a matrix which is optimized for sparse matrices.
 */
class Matrix {
private:
	Graph graph;

	count nRows;
	count nCols;

public:
	/** Default constructor */
	Matrix();

	/**
	 * Constructs the Matrix with size @a dimension x @a dimension.
	 * @param dimension Defines how many rows and columns this matrix has.
	 */
	Matrix(const count &dimension);


	/**
	 * Constructs the Matrix with size @a nRows x @a nCols.
	 * @param nRows Number of rows.
	 * @param nCols Number of columns.
	 */
	Matrix(const count &nRows, const count &nCols);

	/**
	 * Constructs the @a dimension x @a dimension Matrix from the elements at position @a positions with values @values.
	 * @param dimension Defines how many rows and columns this matrix has.
	 * @param positions Defines the position (i,j) of each element specified in @a values.
	 * @param values The values of the matrix elements.
	 */
	Matrix(const count &dimension, const std::vector<std::pair<index, index>> &positions, const std::vector<double> &values);

	/**
	 * Constructs the @a nRows x @a nCols Matrix from the elements at position @a positions with values @values.
	 * @param nRows Defines how many rows this matrix has.
	 * @param nCols Defines how many columns this matrix has.
	 * @param positions Defines the position (i,j) of each element specified in @a values.
	 * @param values The values of the matrix elements.
	 */
	Matrix(const count &nRows, const count &nCols, const std::vector<std::pair<index, index>> &positions, const std::vector<double> &values);

	/**
	 * Constructs the Matrix with the rows in @a rows.
	 * @param rows The rows of the matrix. All rows must have the same dimension.
	 */
	Matrix(const std::vector<Vector> &rows);

	/** Copy constructor */
	Matrix(const Matrix &other);

	/** Default destructor */
	virtual ~Matrix() = default;

	/**
	 * @return Number of rows.
	 */
	count numberOfRows() const;

	/**
	 * @return Number of columns.
	 */
	count numberOfColumns() const;

	/**
	 * @return Value at matrix position (i,j).
	 */
	double operator()(const index &i, const index &j) const;

	/**
	 * Set the matrix at position (@a i, @a j) to @a value.
	 */
	void setValue(const index &i, const index &j, const double &value);

	/**
	 * @return Row @a i of this matrix as vector.
	 */
	Vector row(const index &i) const;

	/**
	 * @return Column @a j of this matrix as vector.
	 */
	Vector column(const index &j) const;

	/**
	 * Adds this matrix to @a other and returns the result.
	 * @return The sum of this matrix and @a other.
	 */
	Matrix operator+(const Matrix &other) const;

	/**
	 * Adds @a other to this matrix.
	 * @return Reference to this matrix.
	 */
	Matrix& operator+=(const Matrix &other);

	/**
	 * Subtracts @a other from this matrix and returns the result.
	 * @return The difference of this matrix and @a other.
	 *
	 */
	Matrix operator-(const Matrix &other) const;

	/**
	 * Subtracts @a other from this matrix.
	 * @return Reference to this matrix.
	 */
	Matrix& operator-=(const Matrix &other);

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar and returns the result.
	 * @return The result of multiplying this matrix with @a scalar.
	 */
	Matrix operator*(const double &scalar) const;

	/**
	 * Multiplies this matrix with a scalar specified in @a scalar.
	 * @return Reference to this matrix.
	 */
	Matrix& operator*=(const double &scalar);

	/**
	 * Multiplies this matrix with @a vector and returns the result.
	 * @return The result of multiplying this matrix with @a vector.
	 */
	Vector operator*(const Vector &vector) const;

	/**
	 * Multiplies this matrix with @a other and returns the result in a new matrix.
	 * @return The result of multiplying this matrix with @a other.
	 */
	Matrix operator*(const Matrix &other) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor and returns the result in a new matrix.
	 * @return The result of dividing this matrix by @a divisor.
	 */
	Matrix operator/(const double &divisor) const;

	/**
	 * Divides this matrix by a divisor specified in @a divisor.
	 * @return Reference to this matrix.
	 */
	Matrix& operator/=(const double &divisor);


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


} /* namespace NetworKit */

template<typename L>
inline void NetworKit::Matrix::forNonZeroElementsInRowOrder(L handle) const {
	for (index i = 0; i < nRows; ++i) {
		graph.forNeighborsOf(i, [&](index j, edgeweight weight){
			handle(i, j, weight);
		});
	}
}

template<typename L>
inline void NetworKit::Matrix::parallelForNonZeroElementsInRowOrder(L handle) const {
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		graph.forNeighborsOf(i, [&](index j, edgeweight weight){
			handle(i, j, weight);
		});
	}
}

template<typename L>
inline void NetworKit::Matrix::parallelForNonZeroElementsInRowOrder(L handle) {
#pragma omp parallel for
	for (index i = 0; i < nRows; ++i) {
		graph.forNeighborsOf(i, [&](index j, edgeweight weight){
			handle(i, j, weight);
		});
	}
}


#endif /* MATRIX_H_ */
