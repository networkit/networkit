/*
 * ParallelMatcher.h
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#ifndef PARALLELMATCHER_H_
#define PARALLELMATCHER_H_

#include "Matcher.h"
#include "../graph/NodeMap.h"

namespace EnsembleClustering {

class ParallelMatcher: public EnsembleClustering::Matcher {

public:

	ParallelMatcher();

	virtual ~ParallelMatcher();

	/**
	 * Apply the parallel matching algorithm described by Manne/Bisseling
	 * Source:  http://link.springer.com/chapter/10.1007%2F978-3-540-68111-3_74?LI=true#page-1
	 *
	 *
	 */
	virtual Matching& run(Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* PARALLELMATCHER_H_ */
