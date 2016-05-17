/*
 * ApproxBetweenness.cpp
 *
 *  Created on: 09.04.2014
 *      Author: cls
 */

#include "ApproxBetweenness.h"
#include "../auxiliary/Random.h"
#include "../distance/Diameter.h"
#include "../graph/Sampling.h"
#include "../graph/Dijkstra.h"
#include "../graph/BFS.h"
#include "../graph/SSSP.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/SignalHandling.h"

#include <math.h>
#include <algorithm>
#include <memory>
#include <omp.h>

namespace NetworKit {

ApproxBetweenness::ApproxBetweenness(const Graph& G, const double epsilon, const double delta, const count diameterSamples, const double universalConstant) : Centrality(G, true), epsilon(epsilon), delta(delta), diameterSamples(diameterSamples), universalConstant(universalConstant) {

}


void ApproxBetweenness::run() {
	Aux::SignalHandler handler;
	scoreData.clear();
	scoreData.resize(G.upperNodeIdBound());

	edgeweight vd = 0;
	if (diameterSamples == 0) {
		INFO("estimating vertex diameter pedantically");
		Diameter diam(G, DiameterAlgo::estimatedPedantic);
		diam.run();
		vd = diam.getDiameter().first;
	} else {
		/**
		* This is an optimization which deviates from the original algorithm.
		* Instead of getting an estimate for each of possibly thousands of connected component and taking the maximum,
		* we sample the graph and take the maximum diameter found. This has a high chance of  hitting the component with the maximum vertex diameter.
		*/
		INFO("estimating vertex diameter roughly");
		Diameter diam(G, DiameterAlgo::estimatedSamples, -1.f, diameterSamples);
		diam.run();
		vd = diam.getDiameter().first;
	}

	INFO("estimated diameter: ", vd);
	r = ceil((universalConstant / (epsilon * epsilon)) * (floor(log2(vd - 2)) + 1 - log(delta)));

	INFO("taking ", r, " path samples");
	// parallelization:
	count maxThreads = omp_get_max_threads();
	DEBUG("max threads: ", maxThreads);
	std::vector<std::vector<double> > scorePerThread(maxThreads, std::vector<double>(G.upperNodeIdBound()));
	DEBUG("score per thread size: ", scorePerThread.size());
	handler.assureRunning();
	#pragma omp parallel for
	for (count i = 1; i <= r; i++) {
		count thread = omp_get_thread_num();
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
		std::unique_ptr<SSSP> sssp;
		if (G.isWeighted()) {
			sssp.reset(new Dijkstra(G, u, true, false, v));
		} else {
			sssp.reset(new BFS(G, u, true, false, v));
		}
		DEBUG("running shortest path algorithm for node ", u);
		if (!handler.isRunning()) continue;
		sssp->run();
		if (!handler.isRunning()) continue;
		if (sssp->numberOfPaths(v) > 0) { // at least one path between {u, v} exists
			DEBUG("updating estimate for path ", u, " <-> ", v);
			// random path sampling and estimation update
			// node s = v;
			node t = v;
			while (t != u)  {
				// sample z in P_u(t) with probability sigma_uz / sigma_us
				std::vector<std::pair<node, double> > choices;
				for (node z : sssp->getPredecessors(t)) {
					bigfloat tmp = sssp->numberOfPaths(z) / sssp->numberOfPaths(t);
					double weight;
					tmp.ToDouble(weight);
					choices.emplace_back(z, weight); 	// sigma_uz / sigma_us
				}
				node z = Aux::Random::weightedChoice(choices);
				assert (z <= G.upperNodeIdBound());
				if (z != u) {
					scorePerThread[thread][z] += 1 / (double) r;
				}
				// s = t;
				t = z;
			}
		}
	}
	handler.assureRunning();

	INFO("adding thread-local scores");
	// add up all thread-local values
	for (auto &local : scorePerThread) {
		G.parallelForNodes([&](node v){
			scoreData[v] += local[v];
		});
	}

	hasRun = true;
}


count ApproxBetweenness::numberOfSamples() {
	INFO("Estimated number of samples", r);
	return r;
}


} /* namespace NetworKit */
