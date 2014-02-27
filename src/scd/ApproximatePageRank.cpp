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

// TODO: pr2 and residual2: Don't create them all the time, use permanent temporary storage
void ApproximatePageRank::push(node u, node seed, std::vector<double>& pr, std::vector<double>& residual)
{
	std::vector<double> pageRank = pr;
	std::vector<double> resid = residual;

	pageRank[u] = pr[u] + alpha * residual[u];
	resid[u] = (1.0 - alpha) * residual[u] / 2.0; // TODO: accelerate by computing constant
	normalizedResid[u] = resid[u] / G.degree(u);

	G.forNeighborsOf(u, [&](node v) {
		resid[v] = residual[v] + (1.0 - alpha) * residual[u] / (2.0 * G.degree(u));  // TODO: accelerate by computing constant
		normalizedResid[v] = resid[v] / G.degree(v);
	});

	pr = pageRank;
	residual = resid;
	TRACE("residual[", u, "]: ", residual[u]);
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

	G.forNodes([&](node v) {
		normalizedResid[v] = residual[v] / (double) G.degree(v);
	});

	auto converged([&](node& argmax) {
		std::vector<double>::iterator max_elem = std::max_element(normalizedResid.begin(), normalizedResid.end());
		argmax = std::distance(normalizedResid.begin(), max_elem);
		TRACE("argmax: ", argmax, ", max: ", (* max_elem), ", eps: ", eps);
		return ((* max_elem) < eps);
	});

	node argMax = seed;
	while (! converged(argMax)) {
		push(argMax, seed, pr, residual);
	}

	return pr;
}

} /* namespace NetworKit */
