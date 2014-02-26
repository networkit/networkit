/*
 * PageRankNibble.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include "PageRankNibble.h"
#include "ApproximatePageRank.h"
#include <cmath>
#include <vector>

namespace NetworKit {

PageRankNibble::PageRankNibble(Graph& g): G(g) {

}

PageRankNibble::~PageRankNibble() {

}

std::set<node> PageRankNibble::bestSweepSet(std::vector<double>& pr) {
	std::set<node> bestCluster;

	// TODO: check conditions, return empty set if not fulfilled for any set


	return bestCluster;
}

std::set<node> PageRankNibble::run(node seed, double phi, unsigned int b) {
	count m = G.numberOfEdges();
	double B = ceil(log2(m));
	double alpha = phi * phi / (225.0 * log(100.0 * sqrt(m)));
	double epsilon = pow(2, -b) / (48.0 * B);

	ApproximatePageRank apr(G, alpha, epsilon);
	std::vector<double> pr = apr.run(seed);

	std::set<node> cluster = bestSweepSet(pr);
	return cluster;
}

} /* namespace NetworKit */
