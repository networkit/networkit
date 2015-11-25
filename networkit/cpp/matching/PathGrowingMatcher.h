/*
 * PathGrowingMatcher.h
 *
 *  Created on: Jun 13, 2013
 *      Author: Henning
 */

#ifndef PATHGROWINGMATCHER_H_
#define PATHGROWINGMATCHER_H_

#include "Matcher.h"
#include "Matching.h"

namespace NetworKit {

/**
 * @ingroup matching
 * Path growing matching algorithm as described by
 *   Hougardy and Drake. TODO: insert DOI
 * Computes an approximate maximum weight matching with guarantee 1/2.
 */
class PathGrowingMatcher: public NetworKit::Matcher {
public:
	/**
	 * @param[in] G Graph for which matching is computed.
	 */
	PathGrowingMatcher(const Graph& G);

	/**
	 * @param[in] G Graph for which matching is computed.
	 */
	PathGrowingMatcher(const Graph& G, const std::vector<double>& edgeScores);

	/**
	 * Runs path growing algorithm to compute approximate maximum weight matching
	 * for graph @a G.
	 * @return Matching (at least half as heavy as maximum weight matching).
	 */
	virtual void run();
};

} /* namespace NetworKit */
#endif /* PATHGROWINGMATCHER_H_ */
