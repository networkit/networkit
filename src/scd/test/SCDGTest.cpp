/*
 * SCDGTest.cpp
 *
 *  Created on: 06.06.2013
 *      Author: cls
 */

#ifndef NOGTEST

#include "SCDGTest.h"

namespace NetworKit {


TEST_F(SCDGTest, testGreedyCommunityExpansion) {
	// TODO: unit test for GreedyCommunityExpansion
}

TEST_F(SCDGTest, testConductance) {

	Graph G(5);
	G.addEdge(0,0);
	G.addEdge(0,1);
	G.addEdge(0,2);
	G.addEdge(0,3);
	G.addEdge(1,2);
	G.addEdge(1,4);
	G.addEdge(2,3);
	G.addEdge(2,4);
	G.addEdge(3,4);

	std::unordered_set<node> first, second, third, fourth, fifth;
	first = {};
	GreedyCommunityExpansion::Conductance conductance1(G, first);
	second = {0};
	GreedyCommunityExpansion::Conductance conductance2(G, second);
	third = {0,1};
	GreedyCommunityExpansion::Conductance conductance3(G, third);
	fourth = {0,1,2};
	GreedyCommunityExpansion::Conductance conductance4(G, fourth);
	fifth = {0,1,2,3};
	GreedyCommunityExpansion::Conductance conductance5(G, fifth);

	double condOne = conductance1.getValue(0);
	double condTwo = conductance2.getValue(1);
	double condThree = conductance3.getValue(2);
	double condFour = conductance4.getValue(3);
	double condFive = conductance5.getValue(4);

	EXPECT_EQ(0.75, 1-condOne) << "1-clustering should have conductance of 0.75";

	EXPECT_GE(0.571429, 1-condTwo) << "2-clustering should have conductance of 4/7";
	EXPECT_LE(0.571428, 1-condTwo) << "2-clustering should have conductance of 4/7";

	EXPECT_GE(0.666667, 1-condThree) << "3-clustering should have conductance of 2/3";
	EXPECT_LE(0.666666, 1-condThree) << "3-clustering should have conductance of 2/3";

	EXPECT_EQ(1, 1-condFour) << "4-clustering should have conductance of 1";

	EXPECT_EQ(1, 1-condFive) << "5-clustering should have conductance of 1";

}

TEST_F(SCDGTest, testBoundarySharpnessTrimming) {

	Graph G (11);
	G.addEdge(0,1);
	G.addEdge(0,2);
	G.addEdge(0,3);
	G.addEdge(0,4);
	G.addEdge(1,2);
	G.addEdge(1,3);
	G.addEdge(1,5);
	G.addEdge(2,4);
	G.addEdge(3,10);
	G.addEdge(4,6);
	G.addEdge(4,7);
	G.addEdge(5,8);
	G.addEdge(5,9);

	std::unordered_set<node> cluster = {0,1,2,3,4,5};
    DummyTrimming trim;
    std::unordered_set<node> community = trim.run(cluster, G);
	EXPECT_EQ(community.size(), 6) << "The community has 6 nodes";

	BoundarySharpness trim2;
	std::unordered_set<node> community2 = trim2.run(cluster, G);
	EXPECT_EQ(community2.size(), 5) << "The community has 5 nodes";

}

TEST_F(SCDGTest, testLocalMdularityM) {

	Graph G(5);
	G.addEdge(0,1);
	G.addEdge(0,2);
	G.addEdge(0,3);
	G.addEdge(0,4);
	G.addEdge(1,2);
	G.addEdge(1,3);
	G.addEdge(1,4);
	G.addEdge(2,3);
	G.addEdge(2,4);
	G.addEdge(3,4);
	G.addEdge(4,4);

	std::unordered_set<node> community = {};
	GreedyCommunityExpansion::LocalModularityM mod (G, community);
	EXPECT_EQ(mod.getValue(0), 0) << "The community should have a local modularity of 0";
	EXPECT_EQ(mod.getValue(4), 0.25) << "The community should have a local modularity of 0.25";

	community.insert(0);
	EXPECT_EQ(mod.getValue(0), 0) << "The community should have a local modularity of 0";
	EXPECT_GE(mod.getValue(1), 0.16666) << "The community should have a local modularity of 1/6";
	EXPECT_LE(mod.getValue(1), 0.16667) << "The community should have a local modularity of 1/6";
	EXPECT_GE(mod.getValue(4), 0.33333) << "The community should have a local modularity of 1/3";
	EXPECT_LE(mod.getValue(4), 0.33334) << "The community should have a local modularity of 1/3";

	community.insert(1);
	EXPECT_EQ(mod.getValue(2), 0.5) << "The community should have a local modularity of 0.5";
	EXPECT_GE(mod.getValue(4), 0.66666) << "The community should have a local modularity of 2/3";
	EXPECT_LE(mod.getValue(4), 0.66667) << "The community should have a local modularity of 2/3";

	community.insert(2);
	EXPECT_EQ(mod.getValue(3), 1.5) << "The community should have a local modularity of 1.5";
	EXPECT_EQ(mod.getValue(4), 1.75) << "The community should have a local modularity of 1.75";

	community.insert(3);
	EXPECT_EQ(mod.getValue(4), 11) << "The community should have a local modularity of 11";

}

//Test with dummy acceptance and conductance and dummy trimming
TEST_F(SCDGTest, testRun) {

	Graph G (12);

	G.addEdge(0,1);
	G.addEdge(0,2);
	G.addEdge(0,3);
	G.addEdge(1,2);
	G.addEdge(1,3);
	G.addEdge(2,3);

	GreedyCommunityExpansion GCE;
	std::unordered_set<node> community = GCE.run(G, 0);
	EXPECT_EQ(community.size(), 2) << "The community should have 2 nodes";

	G.addEdge(0,4);
	G.addEdge(1,5);
	G.addEdge(2,6);
	G.addEdge(3,7);

	community = GCE.run(G, 0);
	EXPECT_EQ(community.size(), 4) << "The community should have 4 nodes";

	community = GCE.run(G, 6);
	EXPECT_EQ(community.size(), 4) << "The community should have 4 nodes";

	G.addEdge(6,8);
	G.addEdge(8,9);
	G.addEdge(8,10);
	G.addEdge(8,11);
	G.addEdge(9,10);
	G.addEdge(9,11);
	G.addEdge(10,11);

	community = GCE.run(G, 0);
	EXPECT_EQ(community.size(), 7) << "The community should have 7 nodes";

	community = GCE.run(G, 4);
	EXPECT_EQ(community.size(), 7) << "The community should have 7 nodes";

	community = GCE.run(G, 6);
	EXPECT_EQ(community.size(), 5) << "The community should have 5 nodes";

	community = GCE.run(G, 8);
	EXPECT_EQ(community.size(), 5) << "The community should have 5 nodes";

	community = GCE.run(G, 9);
	EXPECT_EQ(community.size(), 5) << "The community should have 5 nodes";

}
//
//TEST_F(SCDGTest, testNodeClusterSimilarity) {
//
//	Graph G(5);
//	G.addEdge(0,0);
//	G.addEdge(0,1);
//	G.addEdge(0,2);
//	G.addEdge(0,3);
//	G.addEdge(1,2);
//	G.addEdge(1,4);
//	G.addEdge(2,3);
//	G.addEdge(2,4);
//	G.addEdge(3,4);
//
//
//	std::unordered_set<node> first, second, third, fourth, fifth, shel;
//	first = {0,23};
//	shel = {1,2,3};
//
//	GreedyCommunityExpansion::NodeClusterSimilarity ncs1(G, first, shel);
//	std::cout << ncs1.community << std::endl;
//		std::cout << ncs1.community->size() << std::endl;
//		std::cout << ncs1.shell << std::endl;
//		std::cout << ncs1.shell->size() << std::endl;
//
//	second = {0,1};
//	shel = {2,3,4};
//	GreedyCommunityExpansion::NodeClusterSimilarity ncs2(G, second, shel);
//	third = {0,1,2};
//	shel = {3,4};
//	GreedyCommunityExpansion::NodeClusterSimilarity ncs3(G, third, shel);
//	fourth = {0,1,2,3};
//	shel = {4};
//	GreedyCommunityExpansion::NodeClusterSimilarity ncs4(G, fourth, shel);
//	fifth = {0,1,2,3,4};
//	shel = {};
//	GreedyCommunityExpansion::NodeClusterSimilarity ncs5(G, fifth, shel);
//
//	double ncsOne = ncs1.getValue(1);
//	double ncsTwo = ncs2.getValue(1);
//	std::cout << "[BEGIN] reading graph   " << ncsOne << std::endl;
//
//	double ncsThree = ncs3.getValue(2);
//	std::cout << "[BEGIN] reading graph   " << ncsOne << std::endl;
//
//	double ncsFour = ncs4.getValue(3);
//	std::cout << "[BEGIN] reading graph   " << ncsOne << std::endl;
//
//	double ncsFive = ncs5.getValue(4);
//	std::cout << "[BEGIN] reading graph   " << ncsOne << std::endl;



//	EXPECT_EQ(0, ncsOne) << "1-clustering should have similarity of 0";
//
//	EXPECT_GE(0.571429, ncsTwo) << "2-clustering should have similarity of 4/7";
//	EXPECT_LE(0.571428, ncsTwo) << "2-clustering should have similarity of 4/7";
//
//	EXPECT_GE(0.666667, ncsThree) << "3-clustering should have similarity of 2/3";
//	EXPECT_LE(0.666666, ncsThree) << "3-clustering should have similarity of 2/3";
//
//	EXPECT_EQ(1, ncsFour) << "4-clustering should have similarity of 1";
//
//	EXPECT_EQ(0, ncsFive) << "5-clustering should have similarity of 0";

//}

//
TEST_F(SCDGTest, tryCommunitySubgraph) {
	GraphGenerator gen;
	Graph G = gen.makeCompleteGraph(10);

	node s = 0; // seed node

	GreedyCommunityExpansion GCE;
	std::unordered_set<node> community = GCE.run(G, s);

	// get the subgraph of the community
	Graph sub = Subgraph::fromNodes(G, community);

	// write it to file
//	METISGraphWriter writer;
//	writer.write(sub, "output/CommunitySubgraph.graph");

}

SCDGTest::SCDGTest() {
}

SCDGTest::~SCDGTest() {
}

} /* namespace NetworKit */

#endif /*NOGTEST */
