/*
 * PageRank.h
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#ifndef PAGERANK_H_
#define PAGERANK_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Compute PageRank as node centrality measure.
 */
class PageRank: public NetworKit::Centrality {
protected:
	double damp;
	double tol;

public:
	/**
	 * Constructs the PageRank class for the Graph @a G
	 *
	 * @param[in] G Graph to be processed.
	 * @param[in] damp Damping factor of the PageRank algorithm.
	 * @param[in] tol Error tolerance for PageRank iteration.
	 */
	PageRank(const Graph& G, double damp=0.85, double tol = 1e-8);

	virtual void run();
};

} /* namespace NetworKit */
#endif /* PAGERANK_H_ */
