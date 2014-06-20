/*
 * SimmelianBackboneGTest.cpp
 *
 *  Created on: 31.05.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "SimmelianBackboneGTest.h"

#include "../../backbones/SimmelianBackbone.h"
#include "../../backbones/ChibaNishizekiTriangleCounter.h"

namespace NetworKit {

TEST_F(SimmelianBackboneGTest, testOverlapCounting) {
	//Build up a ranked neighborhood graph. Notation: Alter(Rank)
	//A/100: 1(1),2(1),3(1),5(2),6(2),7(3),9(4)
	//B/200: 1(1),2(2),4(2),6(3),7(3),9(4),10(5),11(6)

	RankedNeighbors neighborsA;
	RankedNeighbors neighborsB;
	neighborsA.push_back(RankedEdge(100,1,0,1));
	neighborsA.push_back(RankedEdge(100,2,0,1));
	neighborsA.push_back(RankedEdge(100,3,0,1));
	neighborsA.push_back(RankedEdge(100,5,0,2));
	neighborsA.push_back(RankedEdge(100,6,0,2));
	neighborsA.push_back(RankedEdge(100,7,0,3));
	neighborsA.push_back(RankedEdge(100,9,0,4));

	neighborsB.push_back(RankedEdge(200,1,0,1));
	neighborsB.push_back(RankedEdge(200,2,0,2));
	neighborsB.push_back(RankedEdge(200,4,0,2));
	neighborsB.push_back(RankedEdge(200,6,0,3));
	neighborsB.push_back(RankedEdge(200,7,0,3));
	neighborsB.push_back(RankedEdge(200,9,0,4));
	neighborsB.push_back(RankedEdge(200,10,0,5));
	neighborsB.push_back(RankedEdge(200,11,0,6));

	SimmelianBackbone simmel(Graph(0), 0.5);


	EXPECT_EQ(1, simmel.getOverlap(neighborsA, neighborsB, 1).overlap) << "wrong overlap";
	EXPECT_DOUBLE_EQ((1.0/3.0), simmel.getOverlap(neighborsA, neighborsB, 1).jaccard) << "wrong jaccard index";

	EXPECT_EQ(1, simmel.getOverlap(neighborsB, neighborsA, 1).overlap) << "wrong overlap";
	EXPECT_DOUBLE_EQ((1.0/3.0), simmel.getOverlap(neighborsB, neighborsA, 1).jaccard) << "wrong jaccard index";

	//Maximum prefix jaccard index for k=4 (5/8).
	EXPECT_EQ(5, simmel.getOverlap(neighborsA, neighborsB, 10).overlap) << "wrong overlap";
	EXPECT_DOUBLE_EQ((5.0/8.0), simmel.getOverlap(neighborsA, neighborsB, 10).jaccard) << "wrong jaccard index";

	EXPECT_EQ(5, simmel.getOverlap(neighborsB, neighborsA, 10).overlap) << "wrong overlap";
	EXPECT_DOUBLE_EQ((5.0/8.0), simmel.getOverlap(neighborsB, neighborsA, 10).jaccard) << "wrong jaccard index";

	EXPECT_EQ(2, simmel.getOverlap(neighborsA, neighborsB, 2).overlap) << "wrong overlap";
	EXPECT_DOUBLE_EQ((2.0/6.0), simmel.getOverlap(neighborsA, neighborsB, 2).jaccard) << "wrong jaccard index";

	EXPECT_EQ(2, simmel.getOverlap(neighborsB, neighborsA, 2).overlap) << "wrong overlap";
	EXPECT_DOUBLE_EQ((2.0/6.0), simmel.getOverlap(neighborsB, neighborsA, 2).jaccard) << "wrong jaccard index";
}

TEST_F(SimmelianBackboneGTest, testRankedNeighborhood) {
	Graph g(10);

	g.addEdge(4,5);
	g.addEdge(4,6);
	g.addEdge(5,6);

	g.addEdge(4,8);
	g.addEdge(4,7);
	g.addEdge(7,8);

	g.addEdge(4,9);
	g.addEdge(8,9);

	//Apply triangle counting algorithm
	ChibaNishizekiTriangleCounter counter;
	edgeCountMap triangles = counter.triangleCounts(g);

	//Actual test: ranked neighborhood
	SimmelianBackbone simmel(Graph(0), 0.5);
	std::vector<RankedNeighbors> neighborhood = simmel.getRankedNeighborhood(g, triangles);

	//Neighborhood of 4
	EXPECT_EQ(5, neighborhood[4].size());
	EXPECT_EQ(RankedEdge(4, 8, 2, 1), neighborhood[4][0]);
	EXPECT_EQ(RankedEdge(4, 5, 1, 2), neighborhood[4][1]);
	EXPECT_EQ(RankedEdge(4, 6, 1, 2), neighborhood[4][2]);
	EXPECT_EQ(RankedEdge(4, 7, 1, 2), neighborhood[4][3]);
	EXPECT_EQ(RankedEdge(4, 9, 1, 2), neighborhood[4][4]);

	//Neighborhood of 8
	EXPECT_EQ(3, neighborhood[8].size());
	EXPECT_EQ(RankedEdge(8, 4, 2, 1), neighborhood[8][0]);
	EXPECT_EQ(RankedEdge(8, 7, 1, 2), neighborhood[8][1]);
	EXPECT_EQ(RankedEdge(8, 9, 1, 2), neighborhood[8][2]);
}

TEST_F(SimmelianBackboneGTest, testRankedNeighborhoodSkippedRanks) {
	Graph g(7);

	g.addEdge(0,1);
	g.addEdge(0,2);
	g.addEdge(0,3);
	g.addEdge(0,4);
	g.addEdge(0,5);
	g.addEdge(0,6);

	g.addEdge(1,2);
	g.addEdge(2,3);
	g.addEdge(3,4);
	g.addEdge(4,5);
	g.addEdge(4,6);

	//Apply triangle counting algorithm
	ChibaNishizekiTriangleCounter counter;
	edgeCountMap triangles = counter.triangleCounts(g);

	//Actual test: ranked neighborhood
	SimmelianBackbone simmel(Graph(0), 0.5);
	std::vector<RankedNeighbors> neighborhood = simmel.getRankedNeighborhood(g, triangles);

	//Neighborhood of 0
	EXPECT_EQ(6, neighborhood[0].size());
	EXPECT_EQ(RankedEdge(0, 4, 3, 1), neighborhood[0][0]);
	EXPECT_EQ(RankedEdge(0, 2, 2, 2), neighborhood[0][1]);
	EXPECT_EQ(RankedEdge(0, 3, 2, 2), neighborhood[0][2]);
	EXPECT_EQ(RankedEdge(0, 1, 1, 4), neighborhood[0][3]);
	EXPECT_EQ(RankedEdge(0, 5, 1, 4), neighborhood[0][4]);
	EXPECT_EQ(RankedEdge(0, 6, 1, 4), neighborhood[0][5]);

	//Neighborhood of 4
	EXPECT_EQ(4, neighborhood[4].size());
	EXPECT_EQ(RankedEdge(4, 0, 3, 1), neighborhood[4][0]);
	EXPECT_EQ(RankedEdge(4, 3, 1, 2), neighborhood[4][1]);
	EXPECT_EQ(RankedEdge(4, 5, 1, 2), neighborhood[4][2]);
	EXPECT_EQ(RankedEdge(4, 6, 1, 2), neighborhood[4][3]);
}

TEST_F(SimmelianBackboneGTest, testOverlapFiltering) {
	Graph g(10);

	g.addEdge(0,1);
	g.addEdge(1,2);
	g.addEdge(2,3);
	g.addEdge(3,4);
	g.addEdge(4,0);

	g.addEdge(0,2);
	g.addEdge(1,3);
	g.addEdge(2,4);
	g.addEdge(3,0);
	g.addEdge(4,1);

	g.addEdge(5,6);
	g.addEdge(6,7);
	g.addEdge(7,8);
	g.addEdge(8,9);
	g.addEdge(9,5);

	g.addEdge(5,7);
	g.addEdge(6,8);
	g.addEdge(7,9);
	g.addEdge(8,5);
	g.addEdge(9,6);

	g.addEdge(0,6);
	g.addEdge(0,5);
	g.addEdge(5,1);

	SimmelianBackbone simmel(g, 2, 1);
	//TODO: ranks are "skipped" in reference implementation. Why?
	Graph b = simmel.calculate();

	EXPECT_EQ(20, b.numberOfEdges());

	EXPECT_FALSE(b.hasEdge(0,6));
	EXPECT_FALSE(b.hasEdge(0,5));
	EXPECT_FALSE(b.hasEdge(5,1));

	EXPECT_TRUE(b.hasEdge(0,1));
	EXPECT_TRUE(b.hasEdge(1,2));
	EXPECT_TRUE(b.hasEdge(2,3));
	EXPECT_TRUE(b.hasEdge(3,4));
	EXPECT_TRUE(b.hasEdge(4,0));
}

TEST_F(SimmelianBackboneGTest, testBackboneTrivial) {
	Graph g(5);

	g.addEdge(0,1);
	g.addEdge(0,2);
	g.addEdge(1,2);

	SimmelianBackbone simmel(g, 1, 0);

	Graph b = simmel.calculate();

	EXPECT_EQ(3, b.numberOfEdges()) << "wrong edge count in backbone";
	EXPECT_EQ(5, b.numberOfNodes()) << "wrong node count in backbone";
}


}
/* namespace NetworKit */

#endif /*NOGTEST */
