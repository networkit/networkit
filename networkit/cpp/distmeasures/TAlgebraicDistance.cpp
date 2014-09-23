/*
 * TAlgebraicDistance.cpp
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#include "../auxiliary/Random.h"
#include "../auxiliary/Log.h"

#include "TAlgebraicDistance.h"

namespace NetworKit {

TAlgebraicDistance::TAlgebraicDistance(const Graph& G) : TNodeDistance(G), numSystems(0), numIters(0), omega(0.0), norm(0) {
	// parameters are set in the initialize method
}

void TAlgebraicDistance::initialize(const Parameters& param) {
	numSystems = param.getInt("numSystems");
	numIters = param.getInt("numIters");
	omega = param.getDouble("omega");
	norm = param.getInt("norm");


	INFO("preprocessing for algebraic distances");
	// random init
	randomInit();

	// main loop
	for (index iter = 0; iter < numIters; ++iter) {
		// store previous iteration
		std::vector<std::vector<double> > oldLoads = loads;

		for (index sys = 0; sys < numSystems; ++sys) {
			G.forNodes([&](node u) {
				double val = 0.0;

				// step 1
				G.forNeighborsOf(u, [&](node v, edgeweight weight) {
					val += weight * oldLoads[sys][v];
				});
				val /= G.weightedDegree(u);

				// step 2
				loads[sys][u] = (1 - omega) * oldLoads[sys][u] + omega * val;
			});
		}
	}
}

double TAlgebraicDistance::distance(node u, node v) {
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

void TAlgebraicDistance::randomInit() {
	count n = G.numberOfNodes();

	// allocate space for loads
	loads.resize(numSystems);
	for (index i = 0; i < numSystems; ++i) {
		loads[i].resize(n);
	}

	for (index i = 0; i < numSystems; ++i) {
		G.forNodes([&](node v) {
			loads[i][v] = Aux::Random::real();
		});
	}
}

} /* namespace NetworKit */
