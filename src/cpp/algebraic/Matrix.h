/*
 * Matrix.h
 *
 *  Created on: 13.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "../graph/Graph.h"
#include <cstdint>
#include "../auxiliary/Log.h"
#include "Vector.h"

class Matrix {
private:
	NetworKit::Graph graph;

public:
	Matrix();
	Matrix(const uint64_t &dimension);
	Matrix(const uint64_t &dimension, const std::vector<std::pair<int, int> > &positions, const std::vector<double> &values);
	Matrix(const uint64_t &dimension, const std::vector<std::pair<NetworKit::node, NetworKit::node> > &positions, const std::vector<double> &values);
	Matrix(const std::vector<Vector> &columns);
	Matrix(const Matrix &other);
	virtual ~Matrix();

	/**
	 * @return Number of rows.
	 */
	uint64_t numberOfRows() const;

	/**
	 * @return Number of columns.
	 */
	uint64_t numberOfColumns() const;

	/**
	 * @return Value at matrix position (i,j).
	 */
	double operator()(const uint64_t &i, const uint64_t &j) const;

	/**
	 * @return Row @a i of this matrix as vector.
	 */
	Vector row(const uint64_t &i) const;

	/**
	 * @return Column @a j of this matrix as vector.
	 */
	Vector column(const uint64_t &j) const;

	/**
	 * Adds this matrix to @a other and returns the result.
	 * Note that the dimensions of the matrices have to be the same.
	 * @return The sum of this matrix and @a other.
	 */
	Matrix operator+(const Matrix &other) const;

	/**
	 * Adds @a other to this matrix.
	 * Note that the dimensions of the matrices have to be the same.
	 * @return Reference to this matrix.
	 */
	Matrix& operator+=(const Matrix &other);

	/**
	 * Subtracts @a other from this matrix and returns the result.
	 * Note that the dimensions of the matrices have to be the same.
	 * @return The difference of this matrix and @a other.
	 *
	 */
	Matrix operator-(const Matrix &other) const;

	/**
	 * Subtracts @a other from this matrix.
	 * Note that the dimensions of the matrices have to be the same.
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
	 * Note that the dimensions must match, i.e. numberOfColumns must be equal
	 * to @a vector dimension.
	 * @return The result of multiplying this matrix with @a vector.
	 */
	Vector operator*(const Vector &vector) const;

	/**
	 * Multiplies this matrix with @a other and returns the result.
	 * Note that the dimensions must match, i.e. numberOfColumns must be equal
	 * to numberOfRows of @a other
	 */
	Matrix operator*(const Matrix &other) const;

	/**
	 * Iterate in parallel over all non-zero elements of the matrix in row order and call handler (lambda closure).
	 */
	template<typename L> void parallelForNonZeroElementsInRowOrder(L handle) const;

	/**
	 * Iterate in parallel over all non-zero elements of the matrix in row order and call handler (lambda closure).
	 */
	template<typename L> void parallelForNonZeroElementsInRowOrder(L handle);
};

template<typename L>
inline void Matrix::parallelForNonZeroElementsInRowOrder(L handle) const {
	auto rowHandle = [&](const NetworKit::node &i) {
		auto elementHandle = [&](const NetworKit::node &j, double value) {
			handle(i, j, value);
		};

		graph.forWeightedNeighborsOf(i, elementHandle);
	};

	graph.parallelForNodes(rowHandle);
}

template<typename L>
inline void Matrix::parallelForNonZeroElementsInRowOrder(L handle) {
	auto rowHandle = [&](const NetworKit::node &i) {
		auto elementHandle = [&](const NetworKit::node &j, double value) {
			handle(i, j, value);
		};

		graph.forWeightedNeighborsOf(i, elementHandle);
	};

	graph.parallelForNodes(rowHandle);
}

#endif /* MATRIX_H_ */
