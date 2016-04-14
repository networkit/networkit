/*
 *
 */

#include "PartitionHubDominance.h"
#include "../auxiliary/SignalHandling.h"
#include <atomic>

void NetworKit::PartitionHubDominance::run() {
	hasRun = false;

	Aux::SignalHandler handler;

	std::vector<std::atomic<count> > maxInternalDeg(P.upperBound());
	std::vector<count> clusterSizes(P.upperBound(), 0);

	handler.assureRunning();

	G.balancedParallelForNodes([&](node u) {
		index c = P[u];

		if (c != none) {
			count internalDeg = 0;
			G.forNeighborsOf(u, [&](node v) {
				if (P[v] == c) {
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

			#pragma omp atomic update
			++clusterSizes[c];
		}
	});

	handler.assureRunning();

	count numClusters = 0;
	weightedAverage = 0;
	unweightedAverage = 0;
	maximumValue = std::numeric_limits<double>::lowest();
	minimumValue = std::numeric_limits<double>::max();
	values.clear();
	values.resize(P.upperBound(), 0);

	for (index i = 0; i < P.upperBound(); ++i) {
		if (clusterSizes[i] > 0) {
			++numClusters;

			double dominance = 1;
			if (clusterSizes[i] > 1) {
				dominance = maxInternalDeg[i] * 1.0 / (clusterSizes[i] - 1);
			}

			values[i] = dominance;
			unweightedAverage += dominance;
			weightedAverage = dominance * clusterSizes[i];

			maximumValue = std::max(dominance, maximumValue);
			minimumValue = std::min(dominance, minimumValue);
		}
	}

	handler.assureRunning();

	unweightedAverage /= numClusters;
	weightedAverage /= G.numberOfNodes();
	hasRun = true;
}
