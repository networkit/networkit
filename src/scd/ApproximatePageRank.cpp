/*
 * ApproximatePageRank.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include "ApproximatePageRank.h"

namespace NetworKit {

ApproximatePageRank::ApproximatePageRank(Graph& g, double alpha_, double epsilon):
		G(g), alpha(alpha_), eps(epsilon)
{

}

void ApproximatePageRank::push(node u, node seed, std::vector<double>& pr, std::vector<double>& residual)
{
	std::vector<double> pr2 = pr;
	std::vector<double> residual2 = residual;

	pr2[u] = pr[u] + alpha * residual[u];
	residual2[u] = (1.0 - alpha) * residual[u] / 2.0; // TODO: accelerate by computing constant

	G.forNeighborsOf(u, [&](node v) {
		residual2[seed] = residual[seed] + (1.0 - alpha) * residual[u] / (2.0 * G.degree(u));  // TODO: accelerate by computing constant
	});

	pr = pr2;
	residual = residual2;
}

ApproximatePageRank::~ApproximatePageRank() {

}

std::vector<double> ApproximatePageRank::run(node seed) {
	// initialize vectors: pr = 0, residual = characteristic(seed)
	count n = G.numberOfNodes();
	std::vector<double> pr(n, 0.0);
	std::vector<double> residual = pr;
	residual[seed] = 1.0;

	std::vector<double> normalizedRes;
	G.forNodes([&](node v) {
		normalizedRes[v] = residual[v] / (double) G.degree(v);
	});

	auto converged([&](node& argmax) {
		std::vector<double>::iterator max_elem = std::max_element(normalizedRes.begin(), normalizedRes.end());
		argmax = std::distance(normalizedRes.begin(), max_elem);
		return ((* max_elem) < eps);
	});

	node argMax = 0;
	while (! converged(argMax)) {
		push(argMax, seed, pr, residual);
	}

	return pr;
}

} /* namespace NetworKit */
