/*
 * TAcceptability.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef TACCEPTABILITY_H_
#define TACCEPTABILITY_H_

#include "TGreedyCommunityExpansion.h"

namespace NetworKit {

class TAcceptability {
public:
	/**
	 * @param[in]	G			pointer to current graph
	 * @param[in]	community	pointer to current community
	 * @param[in]	shell		pointer to current shell
	 */
	TAcceptability(const Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell);

	virtual ~TAcceptability();

	/**
	 * Get the acceptability value for a node. Higher values are better.
	 */
	double getValue(node v);

protected:
	const Graph* G;								//!< pointer to current graph
	std::unordered_set<node>* community;	//!< pointer to current community
	std::unordered_set<node>* shell;		//!< pointer to current shell

};


class TDummyAcceptability: public TAcceptability {

public:

	TDummyAcceptability(const Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell);

	virtual ~TDummyAcceptability();

	double getValue(node v);
};

class TNodeClusterSimilarity: public TAcceptability {

public:

	TNodeClusterSimilarity(const Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell);

	virtual ~TNodeClusterSimilarity();

	 double getValue(node v);
};
} /* namespace NetworKit */
#endif /* TACCEPTABILITY_H_ */
