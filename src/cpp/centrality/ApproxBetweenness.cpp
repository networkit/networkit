/*
 * ApproxBetweenness.cpp
 *
 *  Created on: 09.04.2014
 *      Author: cls
 */

#include "ApproxBetweenness.h"
#include "../auxiliary/Random.h"
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

	count vd = Diameter::estimatedVertexDiameter(G);
	double r = (c / (epsilon * epsilon)) * (floor(log(vd - 2))) + log(1 / delta);

	// double r = (c / (epsilon * epsilon)) * (3 + log(1 / delta));

	for (count i = 1; i <= r; ++i) {
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
			node s, t = v;
			while (t != u) {
				// sample z in P_s(t) with probability sigma_uz / sigma_us
				std::vector<std::pair<node, double> > choices;

				for (node z : path) {
					choices.emplace_back(z, sigma[z] / sigma[s]); 	// sigma_uz / sigma_us
				}
				node z = Aux::Random::weightedChoice(choices);
				if (z != u) {
					scoreData[z] = scoreData[z] + 1 / (double) r;
					s = t;
					t = z;
				}
			}
		}
	}

}


} /* namespace NetworKit */
