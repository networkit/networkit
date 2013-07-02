/*
 * AlgebraicDistances.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */

#include "AlgebraicDistances.h"

namespace NetworKit {

AlgebraicDistances::AlgebraicDistances(const Graph& graph) :
		g(graph) {
	numIters = 0;
	numSystems = 0;
}

AlgebraicDistances::~AlgebraicDistances() {

}

void AlgebraicDistances::randomInit() {
	count n = g.numberOfNodes();
	Aux::RandomProbability randGen;

	// allocate space for loads
	loads.resize(numSystems);
	for (index i = 0; i < numSystems; ++i) {
		loads[i].resize(n);
	}

	for (index i = 0; i < numSystems; ++i) {
		g.forNodes([&](node v) {
			loads[i][v] = randGen.randomFloat();
		});
	}
}

void AlgebraicDistances::preprocess(count numberSystems, count numberIterations,
		double omega) {
	// init
	numSystems = numberSystems;
	numIters = numberIterations;

	// random init
	randomInit();

	// main loop
	for (index iter = 0; iter < numIters; ++iter) {
		// store previous iteration
		std::vector<std::vector<double> > oldLoads = loads;

		for (index sys = 0; sys < numSystems; ++sys) {
			g.forNodes([&](node u) {
				double val = 0.0;

				// step 1
				g.forWeightedNeighborsOf(u, [&](node v, edgeweight weight) {
					val += weight * oldLoads[sys][v];
				});
				val /= g.weightedDegree(u);

				// step 2
				loads[sys][u] = (1 - omega) * oldLoads[sys][u] + omega * val;
			});
		}
	}
}

double AlgebraicDistances::algdist(node u, node v, index norm) const {
	double result = 0.0;

	if (norm == MAX_NORM) { // maximum norm
		for (index sys = 0; sys < numSystems; ++sys) {
			double absDiff = fabs(loads[sys][u] - loads[sys][v]);
			if (absDiff > result) {
				result = absDiff;
			}
		}
	}
	else {
		for (index sys = 0; sys < numSystems; ++sys) {
			double absDiff = fabs(loads[sys][u] - loads[sys][v]);
			result += pow(absDiff, norm);
		}
		result = pow(result, 1.0 / (double) norm);
	}

	return result;
}

double AlgebraicDistances::geometricMeanLoad(node u) const {
	double result = 1.0;
	for (index sys = 0; sys < numSystems; ++sys) {
		result *= loads[sys][u];
	}
	return pow(result, 1.0 / (double) numSystems);
}

} /* namespace NetworKit */
