/*
 * ConcurrentGraphGTest.cpp
 *
 *  Created on: 24.10.2013
 *      Author: cls
 */

#include "ConcurrentGraphGTest.h"
#include "../ConcurrentGraph.h"


namespace NetworKit {

TEST_F(ConcurrentGraphGTest, testParallelEdgeInsertion) {
	count n = 100;
	ConcurrentGraph G(n);

	#pragma omp parallel for
	for (node u = 0; u < n; u++) {
		for (node v = 0; v < n; v++) {
			if (u < v) {
				G.addEdge(u, v);
			}
		}
	}

	count expm = ((n - 1) * n) / 2;
	count m = G.numberOfEdges();
	EXPECT_EQ(expm, m);

	count d = G.degree(0);
	EXPECT_EQ((n - 1), d);

}

} /* namespace NetworKit */
