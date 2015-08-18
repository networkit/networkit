/*
 * KatzCentrality.cpp
 *
 *  Created on: 09.01.2015
 *      Author: Henning
 */

#include "KatzCentrality.h"
#include "../auxiliary/NumericTools.h"

namespace NetworKit {

KatzCentrality::KatzCentrality(const Graph& G, double alpha, double beta, double tol):
		Centrality(G, true), alpha(alpha), beta(beta), tol(tol)
{

}

void KatzCentrality::run() {
	count z = G.upperNodeIdBound();
	std::vector<double> values(z, 1.0);
	scoreData = values;
	double length = 0.0;
	double oldLength = 0.0;

	auto converged([&](double val, double other) {
		// compute residual
		return (Aux::NumericTools::equal(val, other, tol));
	});

	do {
		oldLength = length;

		// iterate matrix-vector product
		G.parallelForNodes([&](node u) {
			values[u] = 0.0;
			// note: inconsistency in definition in Newman's book (Ch. 7) regarding directed graphs
			// we follow the verbal description, which requires to sum over the incoming edges
			G.forInEdgesOf(u, [&](node v, edgeweight ew) {
				values[u] += ew * scoreData[v];
			});
			values[u] *= alpha;
			values[u] += beta;
		});

		// normalize values
		length = 0.0;
		length = G.parallelSumForNodes([&](node u) {
			return (values[u] * values[u]);
		});
		length = sqrt(length);
		G.parallelForNodes([&](node u) {
			values[u] /= length;
		});

//		TRACE("length: ", length);
//		TRACE(values);

		scoreData = values;
	} while (! converged(length, oldLength));

	hasRun = true;

//	// check sign and correct if necessary
//	if (scoreData[0] < 0) {
//		G.parallelForNodes([&](node u) {
//			scoreData[u] = fabs(scoreData[u]);
//		});
//	}
}

} /* namespace NetworKit */
