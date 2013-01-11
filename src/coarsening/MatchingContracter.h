/*
 * MatchingContracter.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef MATCHINGCONTRACTER_H_
#define MATCHINGCONTRACTER_H_

namespace EnsembleClustering {

class MatchingContracter : public Contracter {

public:

	MatchingContracter();

	virtual ~MatchingContracter();

	/**
	 * Contracts graph according to a matching.
	 *
	 * @param[in]	G	fine graph
	 * @param[in]	M	matching
	 *
	 * @param[out]		coarse graph
	 */
	virtual Graph run(Graph& G, Matching& M);
};


} /* namespace EnsembleClustering */
#endif /* MATCHINGCONTRACTER_H_ */
