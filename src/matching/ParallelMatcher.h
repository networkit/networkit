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
#include "../aux/Functions.h"

namespace EnsembleClustering {


/**
 * Parallel matching algorithm as described by Manne/Bisseling
	 * Source:  http://link.springer.com/chapter/10.1007%2F978-3-540-68111-3_74?LI=true#page-1
 */
class ParallelMatcher: public EnsembleClustering::Matcher {
private:
	int attrId; ///< attribute ID of matching scores/weights

public:

	ParallelMatcher(int attrId);

	virtual ~ParallelMatcher();

	virtual Matching run(Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* PARALLELMATCHER_H_ */
