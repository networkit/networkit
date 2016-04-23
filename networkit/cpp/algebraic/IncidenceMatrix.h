/*
 * IncidenceMatrix.h
 *
 *  Created on: 21.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef INCIDENCEMATRIX_H_
#define INCIDENCEMATRIX_H_

#include "../graph/Graph.h"
#include "Matrix.h"
#include "Vector.h"

namespace NetworKit {

/**
 * @ingroup algebraic
 * Incidence matrix of a Graph.
 */
class IncidenceMatrix : public Matrix {
	typedef std::pair<node, node> Edge;

public:
	/**
	 * Constructs the IncidenceMatrix of @a graph.
	 */
	IncidenceMatrix(const Graph &graph);
};

} /* namespace NetworKit */

#endif /* INCIDENCEMATRIX_H_ */
