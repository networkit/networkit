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
	 * Runs path growing algorithm to compute approximate maximum weight matching
	 * for graph @a G.
	 * @param[in] G Graph for which matching is computed.
	 *   All nodes must be alive (static graph).
	 * @return Matching (at least half as heavy as maximum weight matching).
	 */
	virtual Matching run(Graph& G);
};

} /* namespace NetworKit */
#endif /* PATHGROWINGMATCHER_H_ */
