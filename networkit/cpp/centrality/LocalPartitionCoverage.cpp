/*
 *
 */

#include "LocalPartitionCoverage.h"

NetworKit::LocalPartitionCoverage::LocalPartitionCoverage(const NetworKit::Graph &G, const NetworKit::Partition &P): Centrality(G, true, false), P(P) {

}

void NetworKit::LocalPartitionCoverage::run() {
	hasRun = false;

	scoreData.clear();
	scoreData.resize(G.upperNodeIdBound());

	G.balancedParallelForNodes([&](node u) {
		G.forNeighborsOf(u, [&](node, node v, edgeweight ew) {
			if (P[u] == P[v]) scoreData[u] += ew;
		});

		scoreData[u] /= G.weightedDegree(u);
	});

	hasRun = true;
}

double NetworKit::LocalPartitionCoverage::maximum() {
	return 1.0;
}

bool NetworKit::LocalPartitionCoverage::isParallel() const {
	return true;
}

std::string NetworKit::LocalPartitionCoverage::toString() const {
	return "Local partition coverage";
}





