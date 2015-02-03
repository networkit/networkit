/*
 * BarabasiAlbertGenerator.cpp
 *
 *  Created on: May 28, 2013
 *      Author: forigem
 */

#include "../auxiliary/Random.h"

#include "BarabasiAlbertGenerator.h"


namespace NetworKit {

BarabasiAlbertGenerator::BarabasiAlbertGenerator() {
}


BarabasiAlbertGenerator::BarabasiAlbertGenerator(count k,
		count nMax, count n0):k(k), nMax(nMax), n0(n0) {
	if (n0 == 0) {
		n0 = k;
	}
}

Graph BarabasiAlbertGenerator::generate() {
	Graph G = initializeGraph();
	assert (G.numberOfNodes() >= k);
	Aux::ProgressMeter progress(nMax, 200);

	for (count i = n0; i < nMax; i++) {
		count degreeSum = G.numberOfEdges() * 2;
		//DEBUG("Random")ProgressMeter
		node u = G.addNode();
		progress.signal(u);
		std::set<node> targets;
		targets.insert(u);
		int j = 0;
		while (targets.size() - 1 < k) {
			uint64_t random = (uint64_t) Aux::Random::integer(degreeSum);
			j++;
			///if (j > k) throw std::runtime_error("Possible infinite loop detected.");
			bool found = false; // break from node iteration when done
			auto notFound = [&](){ return ! found; };

			G.forNodesWhile(notFound, [&](node v) {


				if (random <= G.degree(v)) {
					found = true; // found a node to connect to
					targets.insert(v);
				}
				random -= G.degree(v);
				//if (j >= G.numberOfNodes() && found==false) throw std::runtime_error("Last node, but still nothing happened.");
			});
		}

		targets.erase(u);

		for (node x : targets) {
			G.addEdge(u, x);
		}

	}

	G.shrinkToFit();
	return G;
}

Graph BarabasiAlbertGenerator::initializeGraph() {
	Graph G(0);

	// initialize the graph with n0 connected nodes
	for (count i = 0; i < n0; i++) {
		node v = G.addNode();
		if (i > 0) G.addEdge(v, v - 1);
	}

	return G;
}


} /* namespace NetworKit */

