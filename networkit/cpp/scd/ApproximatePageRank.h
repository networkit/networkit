/*
 * ApproximatePageRank.h
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#ifndef APPROXIMATEPAGERANK_H_
#define APPROXIMATEPAGERANK_H_

#include <vector>
#include <unordered_map>
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Computes an approximate PageRank vector from a given seed.
 */
class ApproximatePageRank {
protected:
	const Graph& G;
	double alpha;
	double oneMinusAlphaOver2;
	double eps;

	std::unordered_map<node, std::pair<double, double>> pr_res;

	void push(node u, std::set<node>& activeNodes);

public:
	/**
	 * @param g Graph for which an APR is computed.
	 * @param alpha Loop probability of random walk.
	 * @param epsilon Error tolerance.
	 */
	ApproximatePageRank(const Graph& g, double alpha, double epsilon = 1e-12);

	/**
	 * @return Approximate PageRank vector from @a seed with parameters
	 *         specified in the constructor.
	 */
	std::vector<std::pair<node, double>> run(node seed);
};

} /* namespace NetworKit */
#endif /* APPROXIMATEPAGERANK_H_ */
