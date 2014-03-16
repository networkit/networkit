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

class Matrix {
private:
	NetworKit::Graph graph;

public:
	Matrix();
	Matrix(const int &dimension);
	Matrix(const int &dimension, const std::vector<std::pair<NetworKit::node, NetworKit::node> > &pos, const std::vector<double> &values);
	Matrix(const Matrix &other);
	virtual ~Matrix();

	/**
	 * @return Number of rows.
	 */
	int numberOfRows() const;

	/**
	 * @return Number of columns.
	 */
	int numberOfColumns() const;

	/**
	 * @return Value at matrix position (i,j).
	 */
	double operator()(const uint64_t &i, const uint64_t &j) const;

	/**
	 * Adds this matrix to @a other and returns the result.
	 * Note that the dimensions of the matrices have to be the same.
	 * @return The sum of this matrix and @a other.
	 */
	Matrix operator+(const Matrix &other) const;

	/**
	 * Adds @a other to this matrix.
	 * @return Reference to this matrix.
	 */
	Matrix& operator+=(const Matrix &other);
};

#endif /* MATRIX_H_ */
