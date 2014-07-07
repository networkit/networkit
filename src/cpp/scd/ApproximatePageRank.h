/*
 * ApproximatePageRank.h
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#ifndef APPROXIMATEPAGERANK_H_
#define APPROXIMATEPAGERANK_H_

#include <vector>
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup scd
 * Computes an approximate PageRank vector from a given seed.
 */
class ApproximatePageRank {
protected:
	const Graph& G;
	double alpha;
	double oneMinusAlphaOver2;
	double eps;

	std::vector<double> pr; // page rank vector
	std::vector<double> residual;
	std::vector<double> normalizedResid;


	void push(node u, node seed, std::set<node>& activeNodes);

public:
	/**
	 * @param g Graph for which an APR is computed.
	 * @param alpha Loop probability of random walk.
	 * @param epsilon Error tolerance.
	 */
	ApproximatePageRank(Graph& g, double alpha, double epsilon = 1e-12);


	/**
	 * @return Approximate PageRank vector from @a seed with parameters
	 *         specified in the constructor.
	 */
	std::vector<double> run(node seed);
};

} /* namespace NetworKit */
#endif /* APPROXIMATEPAGERANK_H_ */
