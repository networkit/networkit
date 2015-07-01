/*
 * MaxClique.h
 *
 *  Created on: 08.12.2014
 *      Author: Henning
 */

#ifndef MAXCLIQUE_H_
#define MAXCLIQUE_H_

#include "../graph/Graph.h"
#include <set>


namespace NetworKit {

/**
 * Exact algorithm for computing the size of the largest clique in a graph.
 * Worst-case running time is exponential, but in practice the algorithm is fairly fast.
 * Reference: Pattabiraman et al., http://arxiv.org/pdf/1411.7460.pdf
 */
class MaxClique {
protected:
	const Graph& G;
	count maxi;

	/**
	 * Subroutine that goes through every relevant clique containing a certain node in
	 * a recursive fashion and computes the size of the largest.
	 */
	void clique(std::set<node>& U, count size);

public:
	/**
	 * Constructor for maximum clique algorithm.
	 * @param[in] G Graph @a G for which algorithm should be run.
	 * @param[in] lb Lower bound for maximum clique size.
	 */
	MaxClique(const Graph& G, count lb=0);

	/**
	 * Actual maximum clique algorithm. Determines largest clique each vertex
	 * is contained in and returns size of largest. Pruning steps keep running time
	 * acceptable in practice.
	 * @return Size of maximum clique.
	 */
	void run();

	/**
	 * Get size of maximum clique
	 * @return size of maximum clique
	 */
	count getMaxCliqueSize();
};

} /* namespace NetworKit */
#endif /* MAXCLIQUE_H_ */
