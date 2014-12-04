/*
 * ApproximatePageRank.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */


#include <set>
#include <utility>
#include "ApproximatePageRank.h"

namespace NetworKit {

ApproximatePageRank::ApproximatePageRank(const Graph& g, double alpha_, double epsilon):
		G(g), alpha(alpha_), oneMinusAlphaOver2((1.0 - alpha) * 0.5), eps(epsilon)
{

}

void ApproximatePageRank::push(node u, std::set<node>& activeNodes)
{
	double res = pr_res[u].second;
	double mass = oneMinusAlphaOver2 * res / G.degree(u);

	G.forNeighborsOf(u, [&](node v) {
		if (pr_res.find(v) == pr_res.end()) {
			pr_res[v] = std::pair<double, double>(0.0, 0.0);
		}
		pr_res[v] = std::pair<double, double>(pr_res[v].first, pr_res[v].second + mass);
		if ((pr_res[v].second / G.degree(v)) >= eps) {
			activeNodes.insert(v);
		}
	});

	pr_res[u] = std::pair<double, double>(pr_res[u].first + alpha * res, oneMinusAlphaOver2 * res);
	if ((pr_res[u].second / G.degree(u)) >= eps) {
		activeNodes.insert(u);
	}
}


std::vector<std::pair<node, double>> ApproximatePageRank::run(node seed) {
	pr_res[seed] = std::pair<double, double>(0.0, 1.0);
	std::set<node> activeNodes;
	activeNodes.insert(seed);

	while (activeNodes.size() > 0) {
		node v = (* activeNodes.begin());
		activeNodes.erase(v);
		TRACE("queue size: ", activeNodes.size());
		push(v, activeNodes);
	}

	std::vector<std::pair<node, double>> pr;

	for (std::unordered_map<node, std::pair<double, double>>::iterator it = pr_res.begin(); it != pr_res.end(); it++) {
		pr.push_back(std::pair<node, double>(it->first, it->second.first));
	}

	return pr;
}

} /* namespace NetworKit */
