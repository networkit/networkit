/*
 * PageRank.cpp
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#include "PageRank.h"

namespace NetworKit {

NetworKit::PageRank::PageRank(const Graph& G, double damp, double tol):
		Centrality(G, true), damp(damp), tol(tol)
{

}

void NetworKit::PageRank::run() {
	count n = G.numberOfNodes();
	count z = G.upperNodeIdBound();
	double oneOverN = 1.0 / (double) n;
	double teleportProb = (1.0 - damp) / (double) n;
	scoreData.resize(z, oneOverN);
	std::vector<double> pr = scoreData;
	bool isConverged = false;

	while (! isConverged) {
		G.balancedParallelForNodes([&](node u) {
			pr[u] = 0.0;
			G.forNeighborsOf(u, [&](node v) {
				pr[u] += scoreData[v] * G.weight(u, v) / (double) G.weightedDegree(v);
			});
			pr[u] *= damp;
			pr[u] += teleportProb;
		});

		auto converged([&]() {
			double diff = G.parallelSumForNodes([&](node u) {
				double d = scoreData[u] - pr[u];
				return d * d;
			});
			return (sqrt(diff) <= tol);
		});

		isConverged = converged();
		scoreData = pr;
	}

	// make sure scoreData sums up to 1
	double sum = G.parallelSumForNodes([&](node u) {
		return scoreData[u];
	});
	G.parallelForNodes([&](node u) {
		scoreData[u] /= sum;
	});
}

} /* namespace NetworKit */
