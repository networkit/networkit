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

namespace NetworKit {

/**
 * @ingroup algebraic
 * Incidence matrix of a Graph.
 */
class IncidenceMatrix {
	typedef std::pair<node, node> Edge;

private:
	Graph &graph;
	double value(const node &nd, const Edge &edge) const;

public:
	/**
	 * Constructs the IncidenceMatrix of @a graph.
	 */
	IncidenceMatrix(Graph &graph);

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
	 * @return Row @a i of this matrix as vector.
	 */
	Vector row(const index &i) const;

	/**
	 * @return Column @a j of this matrix as vector.
	 */
	Vector column(const index &j) const;

	/**
	 * Multiplies this matrix with @a vector and returns the result.
	 * @return The result of multiplying this matrix with @a vector.
	 */
	Vector operator*(const Vector &vector) const;
};

} /* namespace NetworKit */

#endif /* INCIDENCEMATRIX_H_ */
