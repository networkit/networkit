/*
 * ErdosRenyiGenerator.cpp
 *
 *  Created on: 21.01.2014
 *      Author: Henning
 */

#include "ErdosRenyiGenerator.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/SignalHandling.h"

namespace NetworKit {

ErdosRenyiGenerator::ErdosRenyiGenerator(count nNodes, double prob, bool directed): n(nNodes), p(prob), directed(directed) {

}


/**
 * Returns number of steps you need to wait until the next success (edge) occurs.
 */
static inline count get_next_edge_distance(const double log_cp) {
	return (count) 1 + floor(log(1.0 - Aux::Random::probability()) / log_cp);
}

Graph ErdosRenyiGenerator::generate() {
	Aux::SignalHandler handler;
	Graph G(n, false, directed);
	const double log_cp = log(1.0 - p); // log of counter probability

	// create edges
	node curr = 1;
	node next = -1; // according to Batagelj/Brandes
	while (curr < n) {
		handler.assureRunning();
		// compute new step length
		next += get_next_edge_distance(log_cp);

		// check if at end of row
		while ((next >= curr) && (curr < n)) {
			// adapt to next row
			next = next - curr;
			curr++;
		}

		// insert edge
		if (curr < n) {
			G.addEdge(curr, next);
		}
	}

	G.shrinkToFit();
	return G;
}

} /* namespace NetworKit */
