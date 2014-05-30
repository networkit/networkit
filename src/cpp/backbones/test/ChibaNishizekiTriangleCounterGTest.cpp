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
	EXPECT_EQ(1, (counts[uEdge(0,1)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(0,2)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(1,2)])) << "wrong triangle count";

	EXPECT_EQ(1, (counts[uEdge(0,3)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(3,4)])) << "wrong triangle count";

	EXPECT_EQ(2, (counts[uEdge(0,4)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(0,5)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(4,5)])) << "wrong triangle count";

	//TODO: Remove the following. This is actually a test of uEdge...
	EXPECT_EQ(1, (counts[uEdge(1,0)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(2,0)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(2,1)])) << "wrong triangle count";

	EXPECT_EQ(1, (counts[uEdge(3,0)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(4,3)])) << "wrong triangle count";

	EXPECT_EQ(2, (counts[uEdge(4,0)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(5,0)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(5,4)])) << "wrong triangle count";
}

}
/* namespace NetworKit */

#endif /*NOGTEST */
