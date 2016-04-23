/*
 * MatchingCoarsening.h
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
 * Coarsens graph according to a matching.
 */
class MatchingCoarsening: public GraphCoarsening {

public:
	MatchingCoarsening(const Graph& G, const Matching& M, bool noSelfLoops = false);

	/**
	 * Contracts graph according to a matching.
	 *
	 * @param[in]	G	fine graph
	 * @param[in]	M	matching
	 * @param[in]	noSelfLoops  if true, self-loops are not produced
	 *
	 * @return		coarse graph
	 */
	virtual void run();

private:
	const Matching& M;
	bool noSelfLoops;
};

} /* namespace NetworKit */
#endif /* MATCHINGCONTRACTER_H_ */
