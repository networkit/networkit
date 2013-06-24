/*
 * SCDGTest.cpp
 *
 *  Created on: 06.06.2013
 *      Author: cls
 */

#ifndef NOGTEST

#include "SCDGTest.h"

namespace NetworKit {

SCDGTest::SCDGTest() {
}

SCDGTest::~SCDGTest() {
}
//
//
//TEST_F(SCDGTest, testGreedyCommunityExpansion) {
//	// TODO: unit test for GreedyCommunityExpansion
//}
//
//TEST_F(SCDGTest, testConductance) {
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
//	std::unordered_set<node> first, second, third, fourth, fifth;
//	first = {};
//	Conductance conductance1(G, first);
//	second = {0};
//	Conductance conductance2(G, second);
//	conductance2.volume = 4;
//	conductance2.nBoundaryEdges = 3;
//	third = {0,1};
//	Conductance conductance3(G, third);
//	conductance3.volume = 7;
//	conductance3.nBoundaryEdges = 4;
//	fourth = {0,1,2};
//	Conductance conductance4(G, fourth);
//	conductance4.volume = 11;
//	conductance4.nBoundaryEdges = 4;
//	fifth = {0,1,2,3};
//	Conductance conductance5(G, fifth);
//	conductance5.volume = 14;
//	conductance5.nBoundaryEdges = 3;
//
//	double condOne = conductance1.getValue(0);
//	double condTwo = conductance2.getValue(1);
//	double condThree = conductance3.getValue(2);
//	double condFour = conductance4.getValue(3);
//	double condFive = conductance5.getValue(4);
//
//	EXPECT_EQ(0.75, 1-condOne) << "1-clustering should have conductance of 0.75";
//
//	EXPECT_GE(0.571429, 1-condTwo) << "2-clustering should have conductance of 4/7";
//	EXPECT_LE(0.571428, 1-condTwo) << "2-clustering should have conductance of 4/7";
//
//	EXPECT_GE(0.666667, 1-condThree) << "3-clustering should have conductance of 2/3";
//	EXPECT_LE(0.666666, 1-condThree) << "3-clustering should have conductance of 2/3";
//	EXPECT_EQ(1, 1-condFour) << "4-clustering should have conductance of 1";
//
//	EXPECT_EQ(1, 1-condFive) << "5-clustering should have conductance of 1";
//
//}
//
//TEST_F(SCDGTest, testBoundarySharpnessTrimming) {
//
//	Graph G (11);
//	G.addEdge(0,1);
//	G.addEdge(0,2);
//	G.addEdge(0,3);
//	G.addEdge(0,4);
//	G.addEdge(1,2);
//	G.addEdge(1,3);
//	G.addEdge(1,5);
//	G.addEdge(2,4);
//	G.addEdge(3,10);
//	G.addEdge(4,6);
//	G.addEdge(4,7);
//	G.addEdge(5,8);
//	G.addEdge(5,9);
//
//	std::unordered_set<node> cluster = {0,1,2,3,4,5};
//    DummyTrimming trim;
//    std::unordered_set<node> community = trim.run(cluster, G);
//	EXPECT_EQ(6, community.size()) << "The community has 6 nodes";
//
//	BoundarySharpness trim2;
//	std::unordered_set<node> community2 = trim2.run(cluster, G);
//	EXPECT_EQ(5, community2.size()) << "The community has 5 nodes";
//
//}
//
//TEST_F(SCDGTest, testLocalMdularityM) {
//
//	Graph G(5);
//	G.addEdge(0,1);
//	G.addEdge(0,2);
//	G.addEdge(0,3);
//	G.addEdge(0,4);
//	G.addEdge(1,2);
//	G.addEdge(1,3);
//	G.addEdge(1,4);
//	G.addEdge(2,3);
//	G.addEdge(2,4);
//	G.addEdge(3,4);
//	G.addEdge(4,4);
//
//	std::unordered_set<node> community = {};
//	LocalModularityM mod (G, community);
//	EXPECT_EQ(0, mod.getValue(0)) << "The community should have a local modularity of 0";
//	EXPECT_EQ(0.25, mod.getValue(4)) << "The community should have a local modularity of 0.25";
//
//	community.insert(0);
//	EXPECT_EQ(0, mod.getValue(0)) << "The community should have a local modularity of 0";
//	EXPECT_GE(0.16667, mod.getValue(1)) << "The community should have a local modularity of 1/6";
//	EXPECT_LE(0.16666, mod.getValue(1)) << "The community should have a local modularity of 1/6";
//	EXPECT_GE(0.33334, mod.getValue(4)) << "The community should have a local modularity of 1/3";
//	EXPECT_LE(0.33333, mod.getValue(4)) << "The community should have a local modularity of 1/3";
//
//	community.insert(1);
//	EXPECT_EQ(0.5, mod.getValue(2)) << "The community should have a local modularity of 0.5";
//	EXPECT_GE(0.66667, mod.getValue(4)) << "The community should have a local modularity of 2/3";
//	EXPECT_LE(0.66666, mod.getValue(4)) << "The community should have a local modularity of 2/3";
//
//	community.insert(2);
//	EXPECT_EQ(1.5, mod.getValue(3)) << "The community should have a local modularity of 1.5";
//	EXPECT_EQ(1.75, mod.getValue(4)) << "The community should have a local modularity of 1.75";
//
//	community.insert(3);
//	EXPECT_EQ(11, mod.getValue(4)) << "The community should have a local modularity of 11";
//
//}
//
////Test with dummy acceptance and conductance and dummy trimming
//TEST_F(SCDGTest, testRun) {
//
//	Graph G (12);
//
//	// add clique
//	G.addEdge(0,1);
//	G.addEdge(0,2);
//	G.addEdge(0,3);
//	G.addEdge(1,2);
//	G.addEdge(1,3);
//	G.addEdge(2,3);
//
//	GreedyCommunityExpansion GCE(G);
//	std::unordered_set<node> community= GCE.expandSeed(0);
//	EXPECT_EQ(2, community.size()) << "The community should have 2 nodes";
//
//
//	//add satelites
//	G.addEdge(0,4);
//	G.addEdge(1,5);
//	G.addEdge(2,6);
//	G.addEdge(3,7);
//
//	community = GCE.expandSeed(0);
//
//	EXPECT_EQ(4, community.size()) << "The community should have 4 nodes";
//
//	community = GCE.expandSeed(6);
//	EXPECT_EQ(4, community.size()) << "The community should have 4 nodes";
//
//	// add another clique
//	G.addEdge(6,8);
//	G.addEdge(8,9);
//	G.addEdge(8,10);
//	G.addEdge(8,11);
//	G.addEdge(9,10);
//	G.addEdge(9,11);
//	G.addEdge(10,11);
//
//
//	community = GCE.expandSeed(0);
//	EXPECT_EQ(7, community.size()) << "The community should have 7 nodes";
//
//	community = GCE.expandSeed(4);
//	EXPECT_EQ(7, community.size()) << "The community should have 7 nodes";
//
//	community = GCE.expandSeed(6);
//	EXPECT_EQ(8, community.size()) << "The community should have 8 nodes";
//
//	community = GCE.expandSeed(8);
//	EXPECT_EQ(5, community.size()) << "The community should have 5 nodes";
//
//
//	community = GCE.expandSeed(9);
//	EXPECT_EQ(5, community.size()) << "The community should have 5 nodes";
//
//}
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
//	std::unordered_set<node> first = {0};
//	std::unordered_set<node> shell_1 = {1,2,3};
//	NodeClusterSimilarity ncs1(G, first, shell_1);
//	double ncsOne = ncs1.getValue(1);
//	EXPECT_EQ(0.6, ncsOne) << "1-clustering should have similarity of 0.6";
//
//	std::unordered_set<node> second = {1};
//	std::unordered_set<node> shell_2 = {2,0,4};
//	NodeClusterSimilarity ncs2(G, second, shell_2);
//	double ncsTwo = ncs2.getValue(0);
//	EXPECT_EQ(0.6, ncsTwo) << "2-clustering should have similarity of 0.6";
//
//	std::unordered_set<node> third = {0,1,2};
//	std::unordered_set<node> shell_3 = {3,4};
//	NodeClusterSimilarity ncs3(G, third, shell_3);
//	double ncsThree = ncs3.getValue(3);
//	EXPECT_GE(0.8, ncsThree) << "3-clustering should have similarity of 0.8";
//
//	std::unordered_set<node> fourth = {0,1,2,3};
//	std::unordered_set<node> shell_4 = {4};
//	NodeClusterSimilarity ncs4(G, fourth, shell_4);
//	double ncsFour = ncs4.getValue(4);
//	EXPECT_EQ(0.8, ncsFour) << "4-clustering should have similarity of 0.8";
//
//
//
//	std::unordered_set<node> fifth = {0,1};
//	std::unordered_set<node> shell_5 = {2,3,4};
//	NodeClusterSimilarity ncs5(G, fifth, shell_5);
//	double ncsFive = ncs5.getValue(2);
//	EXPECT_EQ(1, ncsFive) << "5-clustering should have similarity of 1";
//
//}
//
//
//TEST_F(SCDGTest, tryCommunitySubgraph) {
//	METISGraphReader reader;
//	Graph G = reader.read("input/lesmis.graph");
//
//	node s = 0; // seed node
//
//	GreedyCommunityExpansion GCE(G);
//	std::unordered_set<node> community = GCE.expandSeed(s);
//
//	// get the subgraph of the community
//	Graph sub = Subgraph::fromNodes(G, community);
//
//	// write it to file
//	METISGraphWriter writer;
//	writer.write(sub, "output/lesmis-comm0.graph");
//
//}
//
//TEST_F(SCDGTest, testRandomSeedSet) {
//	METISGraphReader reader;
//	Graph G = reader.read("input/jazz.graph");
//
//	RandomSeedSet randSeeds(G);
//
//	count k = 42;
//	std::unordered_set<node> S = randSeeds.getSeeds(k);
//
//	EXPECT_EQ(k, S.size());
//
//	DEBUG("seed set is: " << Aux::setToString(S));
//
//}
//
TEST_F(SCDGTest, testRandomWalkSeedSet) {
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");

	RandomWalkSeedSet walk(G, 2);

	count k = 42;
	std::unordered_set<node> S = walk.getSeeds(k);

	EXPECT_EQ(k, S.size());

	DEBUG("seed set is: " << Aux::setToString(S));

}

TEST_F(SCDGTest, tryGreedyWithSeedSets) {

	METISGraphReader reader;
	Graph G = reader.read("input/caidaROuterLevel.graph");
	Parameters param;
	RandomSeedSet randSeeds(G);

	std::unordered_set<node> seeds = randSeeds.getSeeds(100);
	assert (seeds.size() == 100);

	SelectiveSCAN GCE1(G);
	std::unordered_map<node, std::unordered_set<node>> result1 = GCE1.run(seeds);
	for(auto u :result1) {
		std::cout<<u.first<<"-----------"<<u.second.size()<<std::endl;
	}
}

//TEST_F(SCDGTest, testSelectiveSCAN) {
//	Graph G(10);
//	// add clique
//	G.addEdge(0, 1);
//	G.addEdge(0, 2);
//	G.addEdge(0, 3);
//	G.addEdge(0, 4);
//	G.addEdge(1, 2);
//	G.addEdge(1, 3);
//	G.addEdge(1, 4);
//	G.addEdge(2, 3);
//	G.addEdge(2, 4);
//	G.addEdge(3, 4);
//	G.addEdge(4, 5);
//	G.addEdge(4, 6);
//	G.addEdge(4, 7);
//	G.addEdge(4, 8);
//	G.addEdge(5, 6);
//	G.addEdge(5, 7);
//	G.addEdge(5, 8);
//	G.addEdge(6, 7);
//	G.addEdge(6, 8);
//	G.addEdge(7, 8);
//	G.addEdge(0, 9);
//
//	std::unordered_set<node> seeds = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
//	NeighborhoodDistance distMeasure(G);
//	SelectiveSCAN GCE(G, distMeasure);
//	std::unordered_map<node, std::unordered_set<node>> result = GCE.run(seeds);
//
//	for (auto x : result) {
//		std::cout << "-----------" << x.first << "-----------------"
//				<< std::endl;
//
//		for (node u : x.second) {
//			std::cout << u << std::endl;
//		}
//		std::cout << "----------------------------" << std::endl;
//
//	}
//
//}
//
//TEST_F(SCDGTest, testTSelectiveSCAN) {
//	Graph G(10);
//	// add clique
//	G.addEdge(0, 1);
//	G.addEdge(0, 2);
//	G.addEdge(0, 3);
//	G.addEdge(0, 4);
//	G.addEdge(1, 2);
//	G.addEdge(1, 3);
//	G.addEdge(1, 4);
//	G.addEdge(2, 3);
//	G.addEdge(2, 4);
//	G.addEdge(3, 4);
//	G.addEdge(4, 5);
//	G.addEdge(4, 6);
//	G.addEdge(4, 7);
//	G.addEdge(4, 8);
//	G.addEdge(5, 6);
//	G.addEdge(5, 7);
//	G.addEdge(5, 8);
//	G.addEdge(6, 7);
//	G.addEdge(6, 8);
//	G.addEdge(7, 8);
//	G.addEdge(0, 9);
//
//	std::unordered_set<node> seeds = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
//
//	Parameters param;
//	TSelectiveSCAN<TNeighborhoodDistance> GCE(G, param);
//	std::unordered_map<node, std::unordered_set<node>> result = GCE.run(seeds);
//
//	for (auto x : result) {
//		std::cout << "-----------" << x.first << "-----------------" << std::endl;
//
//		for (node u : x.second) {
//			std::cout << u << std::endl;
//		}
//		std::cout << "----------------------------" << std::endl;
//
//	}
//
//}

//
//TEST_F(SCDGTest, localJaccardTest) {
//	node seedNode = 0;
//	std::unordered_set<node> cluster = {0,1,2,3};
//	Clustering groundTruth(7);
//	JaccardIndex index;
//
//	groundTruth.addToCluster(0, 1);
//	groundTruth.addToCluster(1, 2);
//	groundTruth.addToCluster(0, 4);
//	groundTruth.addToCluster(1, 5);
//
//
//	EXPECT_EQ(0, index.localDissimilarity(seedNode, cluster, groundTruth)) << "The Jaccard index has a value of 0";
//
//	groundTruth.addToCluster(0, 0);
//	cluster = {1,2,3};
//
//	EXPECT_EQ(0, index.localDissimilarity(seedNode, cluster, groundTruth)) << "The Jaccard index has a value of 0";
//
//	cluster = {0,1,2,3};
//	EXPECT_EQ(0.4, index.localDissimilarity(seedNode, cluster, groundTruth)) << "The Jaccard index has a value of 0.4";
//}
//
//TEST_F(SCDGTest, JaccardTest) {
//	std::unordered_map<node,std::unordered_set<node>> communities;
//	communities.insert({1,{1,4}});
//	communities.insert({2,{2,5}});
//	Clustering groundTruth(8);
//	JaccardIndex index;
//
//	groundTruth.addToCluster(0, 1);
//	groundTruth.addToCluster(0, 6);
//	groundTruth.addToCluster(1, 2);
//	groundTruth.addToCluster(1, 4);
//	groundTruth.addToCluster(1, 7);
//
//	EXPECT_LE(0.33333, index.getDissimilarity(communities, groundTruth))<< "The Jaccard index has a value of 1/3";
//	EXPECT_GE(0.33334, index.getDissimilarity(communities, groundTruth))<< "The Jaccard index has a value of 1/3";
//}
//
//TEST_F(SCDGTest, localPrecisionTest) {
//	node seedNode = 0;
//	std::unordered_set<node> cluster = {0,1,2};
//	Clustering groundTruth(5);
//	Precision index;
//
//	groundTruth.addToCluster(0, 0);
//	groundTruth.addToCluster(0, 1);
//	groundTruth.addToCluster(0, 3);
//	groundTruth.addToCluster(1, 4);
//
//	EXPECT_LE(0.66666, index.localDissimilarity(seedNode, cluster, groundTruth))<< "The Precision has a value of 2/3";
//	EXPECT_GE(0.66667, index.localDissimilarity(seedNode, cluster, groundTruth))<< "The Precision has a value of 2/3";
//}
//
//TEST_F(SCDGTest, PrecisionTest) {
//	std::unordered_map<node,std::unordered_set<node>> communities;
//		communities.insert({1,{1,4}});
//		communities.insert({2,{2,5}});
//		Clustering groundTruth(8);
//		Precision index;
//
//		groundTruth.addToCluster(0, 1);
//		groundTruth.addToCluster(0, 6);
//		groundTruth.addToCluster(1, 2);
//		groundTruth.addToCluster(1, 4);
//		groundTruth.addToCluster(1, 7);
//
//	EXPECT_EQ(0.5, index.getDissimilarity(communities, groundTruth))<< "The Precision has a value of 0.5";
//}
//
//TEST_F(SCDGTest, localRecallTest) {
//	node seedNode = 0;
//	std::unordered_set<node> cluster = {0,1,2};
//	Clustering groundTruth(5);
//	Recall index;
//
//	groundTruth.addToCluster(0, 0);
//	groundTruth.addToCluster(0, 1);
//	groundTruth.addToCluster(0, 3);
//	groundTruth.addToCluster(1, 4);
//
//	EXPECT_LE(0.66666, index.localDissimilarity(seedNode, cluster, groundTruth))<< "The Recall has a value of 2/3";
//	EXPECT_GE(0.66667, index.localDissimilarity(seedNode, cluster, groundTruth))<< "The Recall has a value of 2/3";
//}
//
//TEST_F(SCDGTest, RecallTest) {
//	std::unordered_map<node,std::unordered_set<node>> communities;
//		communities.insert({1,{1,4}});
//		communities.insert({2,{2,5}});
//		Clustering groundTruth(8);
//		Recall index;
//
//		groundTruth.addToCluster(0, 1);
//		groundTruth.addToCluster(0, 6);
//		groundTruth.addToCluster(1, 2);
//		groundTruth.addToCluster(1, 4);
//		groundTruth.addToCluster(1, 7);
//
//	EXPECT_EQ(0.4, index.getDissimilarity(communities, groundTruth))<< "The Precision has a value of 0.4";
//}
//TEST_F(SCDGTest, benchmarkGreedy) {
//	METISGraphReader reader;
//	Graph G = reader.read("input/pgp.graph");
//
//	RandomSeedSet randSeeds(G);
//
//	GreedyCommunityExpansion GCE(G);
//	count nRuns = 1000;
//	for (count i = 0; i < nRuns; ++i) {
//		std::unordered_set<node> seeds = randSeeds.getSeeds(1);
//		std::unordered_map<node, std::unordered_set<node>> result = GCE.run(seeds);
//	}
//}

} /* namespace NetworKit */

#endif /*NOGTEST */
