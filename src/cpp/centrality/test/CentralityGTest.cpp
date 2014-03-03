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

 EXPECT_NEAR(0.0, bc[0], 0.001);
 EXPECT_NEAR(0.0, bc[1], 0.001);
 EXPECT_NEAR(15.0, bc[2], 0.001);
 EXPECT_NEAR(3.0, bc[3], 0.001);
 EXPECT_NEAR(3.0, bc[4], 0.001);
 EXPECT_NEAR(1.0, bc[5], 0.001);
}


} /* namespace NetworKit */
