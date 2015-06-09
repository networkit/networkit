/*
 * PageRankNibble.h
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#ifndef PAGERANKNIBBLE_H_
#define PAGERANKNIBBLE_H_

#include <set>
#include <unordered_map>
#include "../graph/Graph.h"
#include "SelectiveCommunityDetector.h"


namespace NetworKit {

/**
 * Variant of PageRank-Nibble algorithm due to Andersen, Chung and Lang.
 * Paper: Local Graph Partitioning using PageRank Vectors.
 * URL: http://www.math.ucsd.edu/~fan/wp/localpartition.pdf
 * Simplifications according to D. Gleich's code at URL https://gist.github.com/dgleich/6201856.
 */
class PageRankNibble : public SelectiveCommunityDetector {

protected:

	double alpha;
	double epsilon;

	std::set<node> bestSweepSet(std::vector<std::pair<node, double>>& pr);

public:
	/**
	 * @param Graph @a g for which PageRank-Nibble is to be performed. Is treated as
	 * unweighted in current implementation.
	 * @param alpha Loop probability of random walk; smaller values tend to produce larger
	 *        communities.
	 * @param eps Tolerance threshold for approximation of PageRank vectors.
	 */
	PageRankNibble(Graph& g, double alpha, double epsilon);

	virtual ~PageRankNibble();

	virtual std::map<node, std::set<node> >  run(std::set<unsigned int>& seeds);


		/**
	 * @param seed Seed node for which a community is to be found.

	 * @return Set of nodes that makes up the best community found around node @a seed.
	 *   If target conductance or target size are not fulfilled, an empty set is returned.
	 */
	std::set<node> expandSeed(node seed);
};

} /* namespace NetworKit */
#endif /* PAGERANKNIBBLE_H_ */
