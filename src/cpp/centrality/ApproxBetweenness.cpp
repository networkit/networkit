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
		if (dijkstra.numberOfPaths(v) > 0) { // at least one path between {u, v} exists
			// random path sampling and estimation update
			node s = v;
			node t = v;
			count i = 0;
			while (t != u)  {
				// DEBUG
				++i;
				if (i >= 50) throw std::runtime_error("too many iterations");
				// DEBUG
				TRACE("u, v, s, t: ", u, " ", v, "  ", s, " ", t);
				// sample z in P_s(t) with probability sigma_uz / sigma_us
				std::vector<std::pair<node, double> > choices;

				for (node z : dijkstra.getPredecessors(t)) {
					choices.emplace_back(z, dijkstra.numberOfPaths(z) / (double) dijkstra.numberOfPaths(s)); 	// sigma_uz / sigma_us
				}
				DEBUG("choices and weightes: ", choices);
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
