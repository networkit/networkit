/*
 * IncidenceMatrix.h
 *
 *  Created on: 21.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef INCIDENCEMATRIX_H_
#define INCIDENCEMATRIX_H_

#include "../graph/Graph.h"
#include "Vector.h"

class IncidenceMatrix {
private:
	NetworKit::Graph &graph;

	typedef std::pair<uint64_t, uint64_t> Edge;
	double value(const uint64_t &node, const Edge &edge) const;
public:
	IncidenceMatrix(NetworKit::Graph &graph);
	virtual ~IncidenceMatrix();

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
	 * Multiplies this matrix with @a vector and returns the result.
	 * @return The result of multiplying this matrix with @a vector.
	 */
	Vector operator*(const Vector &vector) const;
};

#endif /* INCIDENCEMATRIX_H_ */
