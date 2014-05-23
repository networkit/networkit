/*
 * CoarseningGTest.cpp
 *
 *  Created on: 23.05.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "ChibaNishizekiTriangleCounterGTest.h"


#include "../../backbones/ChibaNishizekiTriangleCounter.h"

namespace NetworKit {

TEST_F(ChibaNishizekiTriangleCounterGTest, testTriangleCountsSimple) {
	int64_t n = 6;
	Graph G(n);

	G.addEdge(0,1);
	G.addEdge(0,2);
	G.addEdge(1,2);

	G.addEdge(0,4);
	G.addEdge(0,3);
	G.addEdge(3,4);

	G.addEdge(0,5);
	G.addEdge(4,5);

	EXPECT_EQ(8, G.numberOfEdges()) << "wrong edge count";

	ChibaNishizekiTriangleCounter counter;
	edgeCountMap counts = counter.triangleCounts(G);

	EXPECT_EQ(8, counts.size()) << "wrong triangle count map size";
	EXPECT_EQ(1, (counts[std::pair<node,node>(0,1)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[std::pair<node,node>(0,2)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[std::pair<node,node>(1,2)])) << "wrong triangle count";

	EXPECT_EQ(1, (counts[std::pair<node,node>(0,3)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[std::pair<node,node>(3,4)])) << "wrong triangle count";

	EXPECT_EQ(2, (counts[std::pair<node,node>(0,4)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[std::pair<node,node>(0,5)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[std::pair<node,node>(4,5)])) << "wrong triangle count";
}

} /* namespace NetworKit */

#endif /*NOGTEST */
