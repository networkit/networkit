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
 * Exact algorithm for computing the largest clique in a graph.
 * Worst-case running time is exponential, but in practice the algorithm is fairly fast.
 * Reference: Pattabiraman et al., http://arxiv.org/pdf/1411.7460.pdf
 */
class MaxClique {
protected:
	const Graph& G;
	count maxi;

	/**
	 * Subroutine that
	 */
	void clique(std::set<node>& U, count size);

public:
	/**
	 * Constructor for maximum clique algorithm.
	 * @param[in] G Graph @a G for which algorithm should be run.
	 */
	MaxClique(const Graph& G);

	/**
	 * Actual maximum clique algorithm.
	 * @param[in] lb Lower bound for maximum clique size.
	 * @return Size of maximum clique.
	 */
	count run(count lb = 0);
};

} /* namespace NetworKit */
#endif /* MAXCLIQUE_H_ */
