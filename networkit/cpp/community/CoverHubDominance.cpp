/*
 *
 */

#include "CoverHubDominance.h"
#include "../auxiliary/SignalHandling.h"
#include <atomic>

void NetworKit::CoverHubDominance::run() {
	hasRun = false;
	Aux::SignalHandler handler;

	std::vector<std::atomic<count> > maxInternalDeg(C.upperBound());

	handler.assureRunning();

	G.balancedParallelForNodes([&](node u) {
		for (index c : C[u]) {
			count internalDeg = 0;
			G.forNeighborsOf(u, [&](node v) {
				if (C[v].count(c) > 0) {
					internalDeg++;
				}
			});

			// atomic max operation
			// current maximum value
			count cMaxInt = maxInternalDeg[c].load(std::memory_order_relaxed);

			do {
				if (cMaxInt >= internalDeg) break; // skip if current max is already large enough
			// set new max unless current max has been changed in the meantime, if current max has changed load it (so we can compare again)
			} while (!maxInternalDeg[c].compare_exchange_weak(cMaxInt, internalDeg, std::memory_order_release, std::memory_order_relaxed));
		}
	});

	handler.assureRunning();

	std::vector<count> clusterSizes(C.upperBound(), 0);
	count numMemberships = 0;

	G.forNodes([&](node u) {
		for (index c : C[u]) {
			++clusterSizes[c];
		}

		numMemberships += C[u].size();
	});

	handler.assureRunning();

	unweightedAverage = 0;
	weightedAverage = 0;
	minimumValue = std::numeric_limits<double>::max();
	maximumValue = std::numeric_limits<double>::lowest();
	values.clear();
	values.resize(C.upperBound(), 0);

	count numClusters = 0;

	for (index i = 0; i < C.upperBound(); ++i) {
		if (clusterSizes[i] > 0) {
			++numClusters;

			double dominance = 1;
			if (clusterSizes[i] > 1) {
				dominance = maxInternalDeg[i] * 1.0 / (clusterSizes[i] - 1);
			}

			values[i] = dominance;
			minimumValue = std::min(dominance, minimumValue);
			maximumValue = std::max(dominance, maximumValue);
			unweightedAverage += dominance;
			weightedAverage += dominance * clusterSizes[i];
		}
	}

	handler.assureRunning();

	unweightedAverage /= numClusters;
	weightedAverage /= numMemberships;

	hasRun = true;
}
