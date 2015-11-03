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

namespace NetworKit {


/**
 * @ingroup matching
 * LocalMax matching similar to the one described in the EuroPar13 paper
 * by the Sanders group (Birn, Osipov, Sanders, Schulz, Sitchinava)
 */
class LocalMaxMatcher: public NetworKit::Matcher {
public:

	LocalMaxMatcher(const Graph& G);


	virtual void run();

protected:

};

} /* namespace NetworKit */
#endif /* PARALLELMATCHER_H_ */
