/*
* AllSimplePathsGTest.cpp
*
*  Created on: 27.06.2017
*      Author: Eugenio Angriman
*/

#ifndef NOGTEST

#include "AllSimplePathsGTest.h"
#include "../AllSimplePaths.h"
#include <string>
#include "../../auxiliary/Random.h"


namespace NetworKit {

	TEST_F(AllSimplePathsGTest, testAllSimplePaths) {

		// Generating a graph with 10 nodes and a random number of edges.
		srand (time(NULL));
		int n = 8;
		double p  = 0.25;
		Graph G(n, false, true);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				if (j != i && ((double)rand() / (RAND_MAX)) > p) {
					G.addEdge(i, j);
				}
			}
		}

		INFO("Number of edges: ", G.numberOfEdges());

		AllSimplePaths allSimplePaths(G, 0, 1);
		allSimplePaths.run();

		INFO("Number of simple paths: ", allSimplePaths.numberOfSimplePaths());
		allSimplePaths.parallelForAllSimplePaths([&](std::vector<node> p) {
			INFO(p);
		});
	}

} /* namespace NetworKit */

#endif /*NOGTEST */
