/*
 * AlgebraicMatchingCoarsening.h
 *
 *  Created on: Jul 12, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICMATCHINGCOARSENING_H_
#define NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICMATCHINGCOARSENING_H_

#include "../AlgebraicGlobals.h"
#include "../../coarsening/GraphCoarsening.h"
#include "../../matching/Matching.h"

namespace NetworKit {

/**
 * @ingroup algebraic
 * Implements an algebraic version of the MatchingCoarsening algorithm by computing a projection matrix from fine to
 * coarse.
 */
template<class Matrix>
class AlgebraicMatchingCoarsening : public GraphCoarsening {
public:
	/**
	 * Constructs an instance of AlgebraicMatchingCoarsening for the given Graph @a graph and the corresponding
	 * Matching @a matching. If @a noSelfLoops is set to true (false by default), no self-loops are created
	 * during the coarsening.
	 * @param graph
	 * @param matching
	 * @param noSelfLoops
	 */
	AlgebraicMatchingCoarsening(const Graph& graph, const Matching& matching, bool noSelfLoops = false);

	/**
	 * Computes the coarsening for the graph using the given matching.
	 */
	void run() override;

private:
	Matrix A; // adjacency matrix of the graph
	Matrix P; // projection matrix
	bool noSelfLoops;
};

template<class Matrix>
AlgebraicMatchingCoarsening<Matrix>::AlgebraicMatchingCoarsening(const Graph& graph, const Matching& matching, bool noSelfLoops) : GraphCoarsening(graph), A(Matrix::adjacencyMatrix(graph)), noSelfLoops(noSelfLoops) {
	if (G.isDirected()) throw std::runtime_error("Only defined for undirected graphs.");
	nodeMapping.resize(graph.numberOfNodes());
	std::vector<Triplet> triplets(graph.numberOfNodes());

	count numCoarse = graph.numberOfNodes() - matching.size(graph);
	index idx = 0;
	graph.forNodes([&](node u) {
		index mate = matching.mate(u);
		if ((mate == none) || (u < mate)) {
			// vertex u is carried over to the new level
			nodeMapping[u] = idx;
			++idx;
		}
		else {
			// vertex u is not carried over, receives ID of mate
			nodeMapping[u] = nodeMapping[mate];
		}

		triplets[u] = {u, nodeMapping[u], 1};
	});

	P = Matrix(graph.numberOfNodes(), numCoarse, triplets); // dimensions: fine x coarse
}

template<class Matrix>
void AlgebraicMatchingCoarsening<Matrix>::run() {
	Matrix coarseAdj = P.transpose() * A * P; // Matrix::mTmMultiply performs worse due to high sparsity of P (nnz = n)

	Gcoarsened = Graph(coarseAdj.numberOfRows(), true);
	coarseAdj.forNonZeroElementsInRowOrder([&](node u, node v, double weight) {
		if (u == v && !noSelfLoops) {
			Gcoarsened.addEdge(u, v, weight / 2.0);
		} else if (u < v) {
			Gcoarsened.addEdge(u, v, weight);
		}
	});

	hasRun = true;
}

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICMATCHINGCOARSENING_H_ */
