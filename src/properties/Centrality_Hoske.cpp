/*
 * Centrality_Hoske.cpp
 *
 *  Created on: 05.12.2013
 *      Author: dhoske
 */

#include "Centrality_Hoske.h"

namespace NetworKit {

/** Computes the betweenness centrality of the nodes in G. */
std::vector<double> betweennessCentrality_Hoske(const Graph& G) {
	using namespace std;
	count n = G.numberOfNodes();
	vector<double> betweenness(n);

	return betweenness;
}

} /* namespace NetworKit */
