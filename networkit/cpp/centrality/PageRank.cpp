/*
 * PageRank.cpp
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#include "PageRank.h"
#include "../auxiliary/NumericTools.h"
#include "../auxiliary/SignalHandling.h"

namespace NetworKit {

NetworKit::PageRank::PageRank(const Graph& G, double damp, double tol):
		Centrality(G, true), damp(damp), tol(tol)
{

}

void NetworKit::PageRank::run() {
	Aux::SignalHandler handler;
	count n = G.numberOfNodes();
	count z = G.upperNodeIdBound();
	double oneOverN = 1.0 / (double) n;
	double teleportProb = (1.0 - damp) / (double) n;
	scoreData.resize(z, oneOverN);
	std::vector<double> pr = scoreData;
	bool isConverged = false;

	std::vector<double> deg(z, 0.0);
	G.parallelForNodes([&](node u) {
		deg[u] = (double) G.weightedDegree(u);
	});

	while (! isConverged) {
		handler.assureRunning();
		G.balancedParallelForNodes([&](node u) {
			pr[u] = 0.0;
			G.forInEdgesOf(u, [&](node u, node v, edgeweight w) {
				// note: inconsistency in definition in Newman's book (Ch. 7) regarding directed graphs
				// we follow the verbal description, which requires to sum over the incoming edges
				pr[u] += scoreData[v] * w / deg[v];
			});
			pr[u] *= damp;
			pr[u] += teleportProb;
		});

		auto converged([&]() {
			double diff = G.parallelSumForNodes([&](node u) {
				double d = scoreData[u] - pr[u];
				return d * d;
			});
//			TRACE("sqrt(diff): ", sqrt(diff));
			return (sqrt(diff) <= tol);
		});

		isConverged = converged();
		scoreData = pr;
	}
	handler.assureRunning();
	// make sure scoreData sums up to 1
	double sum = G.parallelSumForNodes([&](node u) {
		return scoreData[u];
	});
	assert(! Aux::NumericTools::equal(sum, 0.0, 1e-15));
	G.parallelForNodes([&](node u) {
		scoreData[u] /= sum;
	});

	hasRun = true;
}

double PageRank::maximum() {
	return 1.0;	// upper bound, could be tighter by assuming e.g. a star graph with n nodes
}

} /* namespace NetworKit */
