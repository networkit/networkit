/*
 * PageRankNibble.h
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#ifndef PAGERANKNIBBLE_H_
#define PAGERANKNIBBLE_H_

#include <set>
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Variant of PageRank-Nibble algorithm due to Andersen, Chung and Lang.
 * Paper: Local Graph Partitioning using PageRank Vectors.
 * URL: http://www.math.ucsd.edu/~fan/wp/localpartition.pdf
 * Simplifications according to D. Gleich's code at URL https://gist.github.com/dgleich/6201856.
 */
class PageRankNibble {
protected:
	Graph& G;

	std::set<node> suitableSweepSet(const std::vector<double>& pr);

	/**
	 * @return Number of elements in @a vec unequal to zero.
	 */
	count supportSize(const std::vector<double>& vec) const;


public:
	/**
	 * @param Graph @a g for which PageRank-Nibble is to be performed. Is treated as
	 * unweighted in current implementation.
	 */
	PageRankNibble(Graph& g);
	virtual ~PageRankNibble();

	/**
	 * @param seed Seed node for which a community is to be found.
	 * @param alpha Loop probability of random walk; smaller values tend to produce larger
	 *        communities.
	 * @param eps Tolerance threshold for approximation of PageRank vectors.
	 * @return Set of nodes that makes up the best community found around node @a seed.
	 *   If target conductance or target size are not fulfilled, an empty set is returned.
	 */
	std::set<node> run(node seed, double alpha, double epsilon);
};

} /* namespace NetworKit */
#endif /* PAGERANKNIBBLE_H_ */
