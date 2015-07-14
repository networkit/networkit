/*
 * MatchingContracter.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MATCHINGCONTRACTER_H_
#define MATCHINGCONTRACTER_H_

#include "GraphCoarsening.h"

#include "../matching/Matching.h"

namespace NetworKit {

/**
 * @ingroup coarsening
 * Contracts graph according to a matching.
 */
class MatchingContracter: public GraphCoarsening {

public:
	/**
	 * Contracts graph according to a matching.
	 *
	 * @param[in]	G	fine graph
	 * @param[in]	M	matching
	 * @param[in]	noSelfLoops  if true, self-loops are not produced
	 *
	 * @return		coarse graph
	 */
	virtual std::pair<Graph, std::vector<node> > run(Graph& G, Matching& M,
			bool noSelfLoops = false);
};

} /* namespace NetworKit */
#endif /* MATCHINGCONTRACTER_H_ */
