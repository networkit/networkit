/*
 * Generator.cpp
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#include "GraphGenerator.h"

namespace EnsembleClustering {

GraphGenerator::GraphGenerator() {
	// TODO Auto-generated constructor stub

}

GraphGenerator::~GraphGenerator() {
	// TODO Auto-generated destructor stub
}

Graph& GraphGenerator::makeErdosRenyiGraph(int64_t n, double p) {
	RandomProbability randP;
	Graph* G = new Graph();
	for (node u = 1; u <= n; ++u) {
		for (node v = u; v <= n; ++v) {
			if (randP.generate() <= p) {
				G->insertEdge(u, v);
			}
		}
	}
	return *G;
}
Graph& GraphGenerator::makeCircularGraph(int64_t n) {
	Graph* G = new Graph();
	for (node u = G->firstNode(); u <= n; ++u) {
		G->insertEdge(u, (u + 1) % n);
	}
	return *G;
}

Graph& GraphGenerator::makeCompleteGraph(int64_t n) {
	RandomProbability randP;
	Graph* G = new Graph();
	for (node u = 1; u <= n; ++u) {
		for (node v = u; v <= n; ++v) {
			G->insertEdge(u, v);
		}
	}
	return *G;
}

} /* namespace EnsembleClustering */
