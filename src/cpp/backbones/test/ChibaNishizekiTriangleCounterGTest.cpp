/*
 * ChibaNishizekiCounterGTest.cpp
 *
 *  Created on: 23.05.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "ChibaNishizekiTriangleCounterGTest.h"

#include "../../backbones/ChibaNishizekiTriangleCounter.h"

namespace NetworKit {

TEST_F(ChibaNishizekiTriangleCounterGTest, testTriangleCountsTrivial) {
	Graph g(5);

	g.addEdge(0,1);
	g.addEdge(0,2);
	g.addEdge(1,2);

	ChibaNishizekiTriangleCounter counter;
	edgeCountMap counts = counter.triangleCounts(g);

	EXPECT_EQ(1, (counts[uEdge(0,1)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(0,2)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(1,2)])) << "wrong triangle count";
}

TEST_F(ChibaNishizekiTriangleCounterGTest, testTriangleCountsSimple) {
	int64_t n = 6;
	Graph g(n);

	g.addEdge(0,1);
	g.addEdge(0,2);
	g.addEdge(1,2);

	g.addEdge(0,4);
	g.addEdge(0,3);
	g.addEdge(3,4);

	g.addEdge(0,5);
	g.addEdge(4,5);

	EXPECT_EQ(8, g.numberOfEdges()) << "wrong edge count";

	ChibaNishizekiTriangleCounter counter;
	edgeCountMap counts = counter.triangleCounts(g);

	EXPECT_EQ(6, g.numberOfNodes()) << "undesired side effect";
	EXPECT_EQ(8, g.numberOfEdges()) << "undesired side effect";

	EXPECT_EQ(8, counts.size()) << "wrong triangle count map size";
	EXPECT_EQ(1, (counts[uEdge(0,1)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(0,2)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(1,2)])) << "wrong triangle count";

	EXPECT_EQ(1, (counts[uEdge(0,3)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(3,4)])) << "wrong triangle count";

	EXPECT_EQ(2, (counts[uEdge(0,4)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(0,5)])) << "wrong triangle count";
	EXPECT_EQ(1, (counts[uEdge(4,5)])) << "wrong triangle count";

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
