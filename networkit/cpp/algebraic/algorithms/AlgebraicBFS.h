/*
 * AlgebraicBFS.h
 *
 *  Created on: Jun 7, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICBFS_H_
#define NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICBFS_H_

#include "../../base/Algorithm.h"
#include "../../graph/Graph.h"
#include "../Vector.h"
#include "../GraphBLAS.h"

namespace NetworKit {

/**
 * @ingroup algebraic
 * Implementation of Breadth-First-Search using the GraphBLAS interface.
 */
template<class Matrix>
class AlgebraicBFS : public Algorithm {
public:
	/**
	 * Constructs an instance of the AlgebraicBFS algorithm for the given Graph @a graph and the given @a source node.
	 * @param graph
	 * @param source
	 */
	AlgebraicBFS(const Graph& graph, node source) : At(Matrix::adjacencyMatrix(graph, MinPlusSemiring::zero()).transpose()), source(source) {}

	/**
	 * Runs a bfs using the GraphBLAS interface from the source node.
	 */
	void run();

	/**
	 * Returns the distance from the source to node @a v.
	 * @param v
	 */
	double distance(node v) const {
		assureFinished();
		assert(v <= At.numberOfRows());
		return distances[v];
	}

private:
	Matrix At;
	node source;
	Vector distances;
};

template<class Matrix>
void AlgebraicBFS<Matrix>::run() {
	count n = At.numberOfRows();
	distances = Vector(n, std::numeric_limits<double>::infinity());
	distances[source] = 0;

	Vector oldDist;
	do {
		oldDist = distances;
		GraphBLAS::MxV<MinPlusSemiring>(At, distances, distances);
	} while (oldDist != distances);

	hasRun = true;
}

} /* namespace NetworKit */



#endif /* NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICBFS_H_ */
