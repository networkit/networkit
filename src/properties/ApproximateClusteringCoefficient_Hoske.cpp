/*
 * ApproximateClusteringCoefficient.cpp
 *
 *  Created on: 14.11.2013
 *      Author: dhoske
 */

#include "ApproximateClusteringCoefficient_Hoske.h"
#include <random>
#include <cmath>

namespace NetworKit { namespace ApproximateClusteringCoefficient {

/* TODO: set seed. */
static std::default_random_engine random;

count niters(double variance, double error) {
	return ceil(-log(error) / (variance * variance));
}

double calculate(bool global, const Graph& G, count k) {
	/* Distribution for choosing vertices: v weighted by d(v) * (d(v) - 1). */
	std::vector<count> weights(G.numberOfNodes());
	G.parallelForNodes([&](node u){
		if (global) {
			weights[u] = std::max(G.degree(u) * (G.degree(u) - 1), count(0));
		} else {
			if (G.degree(u) <= 1) {
				weights[u] = 0;
			} else {
				weights[u] = 1;
			}
		}
	});
	auto random_vertex = std::bind(
		std::discrete_distribution<node>(begin(weights), end(weights)),
		random
	);

	/* Choose random paths of length 2 and determine whether they form a triangle. */
	count ntriangles = 0, npaths = 0;
	for (count i = 0; i < k; ++i) {
		node r = random_vertex();

		/* Nodes with degree <= 1 should be chosen with prob. 0 */
		assert(G.degree(r) >= 2);

		node u = G.randomNeighbor(r);
		node v = u;
		while (u == v) {
			v = G.randomNeighbor(r);
		}

		if (G.hasEdge(u, v)) {
			ntriangles++;
		}
		npaths++;
	}

	if (npaths != 0) {
		return double(ntriangles) / npaths;
	} else {
		return 0.0;
	}
}

}} /* namespace NetworKit */
