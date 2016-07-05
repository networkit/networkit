/*
 * DegreeCentrality.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "DegreeCentrality.h"

namespace NetworKit {

DegreeCentrality::DegreeCentrality(const Graph& G, bool normalized, bool outDeg, bool ignoreSelfLoops) : Centrality(G, normalized), outDeg(outDeg), ignoreSelfLoops(ignoreSelfLoops) {
}

void DegreeCentrality::run() {
	scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);

	if (G.isDirected() && !outDeg) {
		G.parallelForNodes([&](node u) {
			scoreData[u] = G.degreeIn(u);
		});
	} else {
		G.parallelForNodes([&](node u) {
			scoreData[u] = G.degree(u);
		});
	}

	if (normalized) {
		count maxDeg = maximum();
		G.parallelForNodes([&](node u) {
			scoreData[u] = scoreData[u] / maxDeg;
		});
	}
	hasRun = true;
}


double DegreeCentrality::maximum() {
	return G.numberOfNodes() - ignoreSelfLoops;
}


} /* namespace NetworKit */
