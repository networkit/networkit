/*
 * GraphMemoryBenchmark.cpp
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include "GraphMemoryBenchmark.h"

#include "../Graph.h"
#include "../DirectedGraph.h"
#include "../../auxiliary/Random.h"

namespace NetworKit {

TEST_F(GraphMemoryBenchmark, comparision) {
	count n = 100000;
	count m = 20 * n;
	Graph G(n);
	DirectedGraph D(n);
	for (int i = 0; i < m; i++) {
		node u = D.randomNode();
		node v = D.randomNode();
		G.addEdge(u, v);
		D.addEdge(u, v);
	}

	INFO("memory used by Graph instance (in KB): ", G.getMemoryUsage() / 1024);
	G.shrinkToFit();
	INFO("after shrinking: ", G.getMemoryUsage() / 1024);

	INFO("memory used by DirectedGraph instance (in KB): ", D.getMemoryUsage() / 1024);
	D.shrinkToFit();
	INFO("after shrinking:: ", D.getMemoryUsage() / 1024);
}

} /* namespace NetworKit */

#endif /*NOGTEST */
