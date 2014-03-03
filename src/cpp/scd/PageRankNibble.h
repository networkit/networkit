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
 * PageRank-Nibble algorithm due to Andersen, Chung and Lang.
 * Paper: Local Graph Partitioning using PageRank Vectors.
 * URL: http://www.math.ucsd.edu/~fan/wp/localpartition.pdf
 */
class PageRankNibble {
protected:
	Graph& G;

	std::set<node> suitableSweepSet(const std::vector<double>& pr, double phi, unsigned int b, unsigned int B);

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
	 * @param phi Target conductance value, for details see paper.
	 * @param b Target size of returned set in terms of volume, for details see paper.
	 * @return Set of nodes that makes up the best community found around node @a seed.
	 *   If target conductance or target size are not fulfilled, an empty set is returned.
	 */
	std::set<node> run(node seed, double phi, unsigned int b);
};

} /* namespace NetworKit */
#endif /* PAGERANKNIBBLE_H_ */
