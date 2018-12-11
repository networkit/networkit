/*
 * AlgebraicBellmanFord.h
 *
 *  Created on: Jun 6, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICBELLMANFORD_H_
#define NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICBELLMANFORD_H_

#include "../../base/Algorithm.h"
#include "../../graph/Graph.h"
#include "../GraphBLAS.h"

#include <iostream>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Implementation of the Bellman-Ford algorithm using the GraphBLAS interface.
 */
template<class Matrix>
class AlgebraicBellmanFord : public Algorithm {
public:
	/**
	 * Construct an instance of the BellmanFord algorithm for the Graph @a graph and the given @a source node.
	 * @param graph
	 * @param source
	 */
	AlgebraicBellmanFord(const Graph& graph, node source) : At(Matrix::adjacencyMatrix(graph, MinPlusSemiring::zero()).transpose()), source(source), negCycle(false) {}

	/** Default destructor */
	~AlgebraicBellmanFord() = default;

	/**
	 * Compute the shortest path from the source to all other nodes.
	 */
	void run();

	/**
	 * Returns the distance from the source node to @a t.
	 * @param  t Target node.
	 * @return The distance from source to target node @a t.
	 */
	inline edgeweight distance(node t) const {
		assert(t < At.numberOfRows());
		assureFinished();
		return distances[t];
	}

	/**
	 * @return True if there is a negative cycle present in the graph, otherwise false.
	 */
	inline bool hasNegativeCycle() const {
		assureFinished();
		return negCycle;
	}

private:
	const Matrix At;
	node source;
	Vector distances;
	bool negCycle;
};

template<class Matrix>
void AlgebraicBellmanFord<Matrix>::run() {
	count n = At.numberOfRows();
	distances = Vector(n, std::numeric_limits<double>::infinity());
	distances[source] = 0;

	for (index k = 1; k < n; ++k) {
		GraphBLAS::MxV<MinPlusSemiring>(At, distances, distances);
	}

	Vector oldDist = distances;
	GraphBLAS::MxV<MinPlusSemiring>(At, distances, distances);
	negCycle = (oldDist != distances);
	hasRun = true;
}


} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICBELLMANFORD_H_ */
