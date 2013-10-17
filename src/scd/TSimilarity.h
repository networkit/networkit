/*
 * TSimilarity.h
 *
 *  Created on: 14.09.2013
 *      Author: cls, Yassine Marrakchi
 */

#ifndef TSIMILARITY_H_
#define TSIMILARITY_H_

#include "TGreedyExpansion.h"
#include <math.h>

namespace NetworKit {

/**
 * Similarity measures quantify how likely a node from the community shell
 * is to improve the community when it is included.
 */
class TSimilarity {
public:
	/**
	 * @param[in]	G		pointer to current graph
	 */
	TSimilarity(const Graph& G);

	virtual ~TSimilarity();

	/**
	 * Get the Similarity value for a node. Higher values are better.
	 *
	 * @param[in]	u 		first node of the considered pair
	 * @param[in]	v		second node of the considered pair
	 */
	double getValue(node u, node v);

protected:
	const Graph* G;								//!< pointer to current graph
};

/**
 * Return the same value for every node.
 *
 * This class is used for variants without Similarity
 */
class TDummySimilarity: public TSimilarity {

public:

	TDummySimilarity(const Graph& G);

	virtual ~TDummySimilarity();

	double getValue(node u, node v);
};

/**
 * Get the node cluster similarity value for a node.
 *
 */
class TNodesSimilarity: public TSimilarity {

public:

	TNodesSimilarity(const Graph& G);

	virtual ~TNodesSimilarity();

	 double getValue(node u, node v);
};
} /* namespace NetworKit */
#endif /* TSIMILARITY_H_ */
