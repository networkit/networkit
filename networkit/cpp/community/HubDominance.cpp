/*
 * HubDominance implementation
 *
 * Created: 2014-08-13
 * Author: Michael Hamann
 */

#include "HubDominance.h"

double NetworKit::HubDominance::getQuality(const NetworKit::Partition &zeta, const NetworKit::Graph &G) {
	std::vector<count> maxInternalDeg(zeta.upperBound());
	auto clusterSizes = zeta.subsetSizeMap();

	G.forNodes([&](node u) {
		index c = zeta[u];

		if (c != none) {
			count internalDeg = 0;
			G.forNeighborsOf(u, [&](node v) {
				if (zeta[v] == c) {
					internalDeg++;
				}
			});
			maxInternalDeg[c] = std::max(maxInternalDeg[c], internalDeg);
		}
	});

	double dominance = 0;

	for (auto it : clusterSizes) {
		if (it.second > 1) {
			dominance += maxInternalDeg[it.first] * 1.0 / (it.second - 1);
		} else {
			dominance += 1;
		}
	}

	return dominance / clusterSizes.size();
}

double NetworKit::HubDominance::getQuality(const NetworKit::Cover &zeta, const NetworKit::Graph &G) {
	std::vector<count> maxInternalDeg(zeta.upperBound());
	auto clusterSizes = zeta.subsetSizeMap();

	G.parallelForNodes([&](node u) {
		for (index c : zeta[u]) {
			count internalDeg = 0;
			G.forNeighborsOf(u, [&](node v) {
				if (zeta[v].count(c) > 0) {
					internalDeg++;
				}
			});

			if (maxInternalDeg[c] < internalDeg) { // assumption: we won't be in this if very often, so it's okay to use pragma omp critical here
				#pragma omp critical
				maxInternalDeg[c] = std::max(maxInternalDeg[c], internalDeg);
			}
		}
	});

	double dominance = 0;

	for (auto it : clusterSizes) {
		if (it.second > 1) {
			dominance += maxInternalDeg[it.first] * 1.0 / (it.second - 1);
		} else {
			dominance += 1;
		}
	}

	return dominance / clusterSizes.size();
}
