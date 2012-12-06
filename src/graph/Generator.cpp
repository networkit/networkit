/*
 * Generator.cpp
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#include "Generator.h"

namespace EnsembleClustering {

Generator::Generator() {
	// TODO Auto-generated constructor stub

}

Generator::~Generator() {
	// TODO Auto-generated destructor stub
}

Graph& Generator::makeErdosRenyiGraph(int64_t n, double p) {
	RandomProbability randP;
	Graph G;
	for (node u = 1; u <= n; ++u) {
		for (node v = u; v <= n; ++v) {
			if (randP.generate() <= p) {
				G.insertEdge(u, v);
			}
		}
	}
	return G;
}

Graph& Generator::makeCircularGraph(int64_t n) {
	Graph G;
	for (node u = G.firstNode(); u <= n; ++u) {
		G.insertEdge(u, (u + 1) % n);
	}
	return G;
}

} /* namespace EnsembleClustering */
