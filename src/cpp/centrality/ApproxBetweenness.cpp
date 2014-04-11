/*
 * ApproxBetweenness.cpp
 *
 *  Created on: 09.04.2014
 *      Author: cls
 */

#include "ApproxBetweenness.h"
#include "../properties/Diameter.h"
#include "../graph/Sampling.h"
#include "../graph/Dijkstra.h"

#include <math.h>

namespace NetworKit {

ApproxBetweenness::ApproxBetweenness(const Graph& G, bool normalized, double epsilon, double delta) : Centrality(G, normalized), epsilon(epsilon), delta(delta) {

}


void ApproxBetweenness::run() {
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);

	double c = 1; // TODO: what is c?

	// TODO: get vertex diameter
	count vd = Diameter::estimatedVertexDiameter(G);

	double r = (c / (epsilon * epsilon)) * (floor(log(vd - 2))) + log(1 / delta);

	for (count i = 0; i <= r; ++i) {
		// sample random node pair
		node u, v;
		u = Sampling::randomNode(G);
		do {
			v = Sampling::randomNode(G);
		} while (v == u);
		Dijkstra dijkstra(G, u);
		dijkstra.run();
		std::vector<edgeweight> sigma = dijkstra.getDistances();
		std::vector<node> path = dijkstra.getPath(v); 	// this selects one shortest path - there may be several
		if (path.size() > 0) { // path exists
			// random path sampling and estimation update
			node j, s, t = v;
			while (t != u) {
				// TODO: sample z in P_s(t) with probability sigma_uz / sigma_us
				if (z != u) {
					// TODO:
				}
			}
		}
	}

}


} /* namespace NetworKit */
