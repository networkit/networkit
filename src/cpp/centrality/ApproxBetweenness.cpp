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
#include "../graph/BFS.h"
#include "../graph/SSSP.h"

#include <math.h>

namespace NetworKit {

ApproxBetweenness::ApproxBetweenness(const Graph& G, double epsilon, double delta) : Centrality(G, true), epsilon(epsilon), delta(delta) {

}


void ApproxBetweenness::run() {
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);

	double c = 1; // TODO: what is c?

	INFO("estimating vertex diameter");
	count vd = Diameter::estimatedVertexDiameter(G);
	double r = (c / (epsilon * epsilon)) * (floor(log(vd - 2))) + log(1 / delta);

	// double r = (c / (epsilon * epsilon)) * (3 + log(1 / delta));

	INFO("taking ", r, " path samples");
	for (count i = 1; i <= r; ++i) {
		DEBUG("sample ", i);
		// if (i >= 1000) throw std::runtime_error("too many iterations");
		// DEBUG
		// sample random node pair
		node u, v;
		u = Sampling::randomNode(G);
		do {
			v = Sampling::randomNode(G);
		} while (v == u);

		// runs faster for unweighted graphs
		SSSP* sssp;
		if (G.isWeighted()) {
			sssp = new Dijkstra(G, u);
		} else {
			sssp = new BFS(G, u);
		}
		DEBUG("running Dijkstra for node ", u);
		sssp->run();
		if (sssp->numberOfPaths(v) > 0) { // at least one path between {u, v} exists
			DEBUG("updating estimate for path ", u, " <-> ", v);
			// random path sampling and estimation update
			node s = v;
			node t = v;
			while (t != u)  {
				// sample z in P_u(t) with probability sigma_uz / sigma_us
				std::vector<std::pair<node, double> > choices;

				for (node z : sssp->getPredecessors(t)) {
					choices.emplace_back(z, sssp->numberOfPaths(z) / (double) sssp->numberOfPaths(s)); 	// sigma_uz / sigma_us
				}
				node z = Aux::Random::weightedChoice(choices);
				if (z != u) {
					scoreData[z] = scoreData[z] + 1 / (double) r;
				}
				s = t;
				t = z;
			}
		}

		delete sssp; // free heap memory
	}

}


} /* namespace NetworKit */
