/*
 * MatchingContracter.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MATCHINGCONTRACTER_H_
#define MATCHINGCONTRACTER_H_

#include "Contracter.h"

#include "../matching/Matching.h"

namespace NetworKit {

/**
 * Contracts graph according to a matching.
 */
class MatchingContracter: public Contracter {

public:

	MatchingContracter();

	virtual ~MatchingContracter();

	/**
	 * Contracts graph according to a matching.
	 *
	 * @param[in]	G	fine graph
	 * @param[in]	M	matching
	 * @param[in]	noSelfLoops  if true, self-loops are not produced
	 *
	 * @return		coarse graph
	 */
	virtual std::pair<Graph, NodeMap<node> > run(Graph& G, Matching& M,
			bool noSelfLoops = false);
};

} /* namespace NetworKit */
#endif /* MATCHINGCONTRACTER_H_ */
