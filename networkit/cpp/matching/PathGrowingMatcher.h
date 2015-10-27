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
 * TODO: class documentation
 */
class PathGrowingMatcher: public NetworKit::Matcher {
protected:

public:

	virtual Matching run(const Graph& G, const bool matchSelfLoops = true);
};

} /* namespace NetworKit */
#endif /* PATHGROWINGMATCHER_H_ */
