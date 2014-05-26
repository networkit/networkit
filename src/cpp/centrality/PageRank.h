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
 * Compute PageRank as node centrality measure.
 */
class PageRank: public NetworKit::Centrality {
protected:
	double damp;
	double tol;

public:
	/**
	 * @param[in] G Graph to be processed.
	 * @param[in] damp Damping factor of the PageRank algorithm.
	 * @param[in] tol Error tolerance for PageRank iteration.
	 */
	PageRank(Graph& G, double damp, double tol = 1e-9);

	virtual void run();
};

} /* namespace NetworKit */
#endif /* PAGERANK_H_ */
