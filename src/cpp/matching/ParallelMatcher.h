/*
 * ParallelMatcher.h
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef PARALLELMATCHER_H_
#define PARALLELMATCHER_H_

#include <set>
#include <algorithm>

#include "Matcher.h"
#include "../graph/NodeMap.h"

namespace NetworKit {


/**
 * LocalMax matching as described in the EuroPar13 paper by the Sanders group
 */
class LocalMaxMatcher: public NetworKit::Matcher {
private:
	uint64_t attrId; ///< attribute ID of matching scores/weights

public:

	LocalMaxMatcher(uint64_t attrId);

	virtual ~LocalMaxMatcher();

	virtual Matching run(Graph& G);
};

} /* namespace NetworKit */
#endif /* PARALLELMATCHER_H_ */
