/*
 * DegreeCentrality.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "DegreeCentrality.h"
#include "../properties/GraphProperties.h"

namespace NetworKit {

DegreeCentrality::DegreeCentrality(const Graph& G, bool normalized) : Centrality(G, normalized) {
}

void DegreeCentrality::run() {
	scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);

	G.parallelForNodes([&](node u) {
		scoreData[u] = G.degree(u);
	});

	if (normalized) {
		count maxDeg = GraphProperties::minMaxDegree(G).second;
		G.parallelForNodes([&](node u) {
			scoreData[u] = scoreData[u] / maxDeg;
		});
	}

	hasRun = true;
}


double DegreeCentrality::maximum() {
	return G.numberOfNodes();
}


} /* namespace NetworKit */
