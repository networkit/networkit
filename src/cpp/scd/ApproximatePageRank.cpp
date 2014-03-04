/*
 * ApproximatePageRank.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include "ApproximatePageRank.h"

namespace NetworKit {

ApproximatePageRank::ApproximatePageRank(Graph& g, double alpha_, double epsilon):
		G(g), alpha(alpha_), oneMinusAlphaOver2((1.0 - alpha) * 0.5), eps(epsilon)
{
	count n = G.numberOfNodes();
	pageRank.resize(n);
	resid.resize(n);
}

void ApproximatePageRank::push(node u, node seed, std::vector<double>& pr, std::vector<double>& residual)
{
	pageRank[u] = pr[u] + alpha * residual[u];
	resid[u] = oneMinusAlphaOver2 * residual[u];
	normalizedResid[u] = resid[u] / G.degree(u);
//	TRACE("normalizedResid[", u, "]: ", normalizedResid[u]);
	double mass = oneMinusAlphaOver2 * residual[u] / G.degree(u);

	G.forNeighborsOf(u, [&](node v) {
		resid[v] = residual[v] + mass;
		normalizedResid[v] = resid[v] / G.degree(v);
//		TRACE("normalizedResid[", v, "]: ", normalizedResid[v]);
	});

	pr = pageRank;
	residual = resid;
//	TRACE("residual[", u, "]: ", residual[u]);
}

ApproximatePageRank::~ApproximatePageRank() {

}

std::vector<double> ApproximatePageRank::run(node seed) {
	// initialize vectors: pr = 0, residual = characteristic(seed)
	count n = G.numberOfNodes();
	std::vector<double> pr(n, 0.0);
	resid = pr;
	normalizedResid = pr;
	std::vector<double> residual = pr;
	residual[seed] = 1.0;
	normalizedResid[seed] = 1.0 / (double) G.degree(seed);


	auto converged([&](node& argmax) {
		// TODO: accelerate
		std::vector<double>::iterator max_elem = std::max_element(normalizedResid.begin(), normalizedResid.end());
		argmax = std::distance(normalizedResid.begin(), max_elem);
//		TRACE("argmax: ", argmax, ", max: ", (* max_elem), ", eps: ", eps);
		return ((* max_elem) < eps);
	});


	node argMax = seed;
	while (! converged(argMax)) {
		push(argMax, seed, pr, residual);
	}

	return pr;
}

} /* namespace NetworKit */
