/*
 * CentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "CentralityGTest.h"
#include "../Betweenness.h"

namespace NetworKit {

TEST_F(CentralityGTest, testBetweenness) {
 /* Graph:
    0    3
     \  / \
      2    5
     /  \ /
    1    4
 */
	count n = 6;
	Graph G(n);

	G.addEdge(0, 2);
	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 5);
	G.addEdge(4, 5);

	Betweenness centrality = Betweenness(G);
	centrality.run();
	std::vector<double> bc = centrality.scores();

	const double tol = 1e-3;
	EXPECT_NEAR(0.0, bc[0], tol);
	EXPECT_NEAR(0.0, bc[1], tol);
	EXPECT_NEAR(15.0, bc[2], tol);
	EXPECT_NEAR(3.0, bc[3], tol);
	EXPECT_NEAR(3.0, bc[4], tol);
	EXPECT_NEAR(1.0, bc[5], tol);
}

TEST_F(CentralityGTest, testBetweennessWeighted) {
 /* Graph:
    0    3   6
     \  / \ /
      2 -- 5
     /  \ / \
    1    4   7

    Edges in the upper row have weight 3,
    the edge in the middle row has weight 1.5,
    edges in the lower row have weight 2.
 */
	count n = 8;
	Graph G(n, true);

	G.addEdge(0, 2, 3);
	G.addEdge(1, 2, 2);
	G.addEdge(2, 3, 3);
	G.addEdge(2, 4, 2);
	G.addEdge(2, 5, 1.5);
	G.addEdge(3, 5, 3);
	G.addEdge(4, 5, 2);
	G.addEdge(5, 6, 3);
	G.addEdge(5, 7, 2);

	Betweenness centrality = Betweenness(G);
	centrality.run();
	std::vector<double> bc = centrality.scores();

	const double tol = 1e-3;
	EXPECT_NEAR(0.0, bc[0], tol);
	EXPECT_NEAR(0.0, bc[1], tol);
	EXPECT_NEAR(23.0, bc[2], tol);
	EXPECT_NEAR(0.0, bc[3], tol);
	EXPECT_NEAR(0.0, bc[4], tol);
	EXPECT_NEAR(23.0, bc[5], tol);
	EXPECT_NEAR(0.0, bc[6], tol);
	EXPECT_NEAR(0.0, bc[7], tol);
}


} /* namespace NetworKit */
