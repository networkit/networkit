/*
 * ApproximatePageRank.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include <set>

#include "ApproximatePageRank.h"

namespace NetworKit {

ApproximatePageRank::ApproximatePageRank(Graph& g, double alpha_, double epsilon):
		G(g), alpha(alpha_), oneMinusAlphaOver2((1.0 - alpha) * 0.5), eps(epsilon)
{

}

void ApproximatePageRank::push(node u, node seed, std::set<node>& activeNodes)
{
	double res = residual[u];
	pr[u] = pr[u] + alpha * res;
//	TRACE("normalizedResid[", u, "]: ", normalizedResid[u]);
	double mass = oneMinusAlphaOver2 * res / G.degree(u);

	G.forNeighborsOf(u, [&](node v) {
		residual[v] = residual[v] + mass;
		normalizedResid[v] = residual[v] / G.degree(v);
		if (normalizedResid[v] >= eps) {
			activeNodes.insert(v);
		}
//		TRACE("normalizedResid[", v, "]: ", normalizedResid[v]);
	});

	residual[u] = oneMinusAlphaOver2 * res;
	normalizedResid[u] = residual[u] / G.degree(u);
	if (normalizedResid[u] >= eps) {
		activeNodes.insert(u);
	}

//	TRACE("normalizedResidual[", u, "]: ", normalizedResid[u]);
}

std::vector<double> ApproximatePageRank::run(node seed) {
	count n = G.upperNodeIdBound();
	pr.assign(n, 0.0);
	residual = pr;
	residual[seed] = 1.0;
	normalizedResid = pr;
	normalizedResid[seed] = 1.0 / (double) G.degree(seed);

	std::set<node> activeNodes;
	activeNodes.insert(seed);

	while (activeNodes.size() > 0) {
		node v = (* activeNodes.begin());
		activeNodes.erase(v);
//		TRACE("queue size: ", activeNodes.size());
		push(v, seed, activeNodes);
	}

	return pr;
}

} /* namespace NetworKit */
