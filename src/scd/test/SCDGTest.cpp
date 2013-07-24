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

TEST_F(SCDGTest, tryCommunitySubgraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/lesmis.graph");

	node s = 0; // seed node
	std::unordered_set<node> tmp = {};
	std::unordered_map<node, count> bound =  {};
	DummySimilarity similarity(G, tmp, tmp);
	Conductance objective(G, tmp, bound);
	DummyTrimming trimming;

	GreedyCommunityExpansion GCE(G, similarity, objective, trimming);

	std::unordered_set<node> community = GCE.expandSeed(s);

	// get the subgraph of the community
	Graph sub = Subgraph::fromNodes(G, community);

	// write it to file
	METISGraphWriter writer;
	writer.write(sub, "output/lesmis-comm0.graph");

}

TEST_F(SCDGTest, testRandomSeedSet) {
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");

	RandomSeedSet randSeeds(G);

	count k = 42;
	std::unordered_set<node> S = randSeeds.getSeeds(k);

	EXPECT_EQ(k, S.size());

	DEBUG("seed set is: " << Aux::setToString(S));

}

TEST_F(SCDGTest, testRandomWalkSeedSet) {
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");

	RandomWalkSeedSet walk(G, 2);

	count k = 42;
	std::unordered_set<node> S = walk.getSeeds(k);

	EXPECT_EQ(k, S.size());

	DEBUG("seed set is: " << Aux::setToString(S));

}

TEST_F(SCDGTest, benchmarkGreedy) {
	METISGraphReader reader;
	Graph G = reader.read("input/pgp.graph");

	RandomSeedSet randSeeds(G);

	std::unordered_set<node> tmp = {};
	std::unordered_map<node, count> bound = {};
	DummySimilarity similarity(G, tmp, tmp);
	LocalModularityM objective(G, tmp, bound);
	BoundarySharpness trimming;

	GreedyCommunityExpansion GCE(G, similarity, objective, trimming);
	count nRuns = 1000;
	for (count i = 0; i < nRuns; ++i) {
		std::unordered_set<node> seeds = randSeeds.getSeeds(5);
		std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> result = GCE.run(seeds);
	}
}

// Test Dissimilarity Measures
TEST_F(SCDGTest, localJaccardTest) {
	node seedNode = 0;
	std::unordered_set<node> cluster = { 0, 1, 2, 3 };
	Clustering groundTruth(7);
	JaccardIndex index;

	groundTruth.addToCluster(0, 1);
	groundTruth.addToCluster(1, 2);
	groundTruth.addToCluster(0, 4);
	groundTruth.addToCluster(1, 5);

	EXPECT_EQ(0, index.localDissimilarity(seedNode, cluster, groundTruth))
			<< "The Jaccard index has a value of 0";

	groundTruth.addToCluster(0, 0);
	cluster = {1,2,3};

	EXPECT_EQ(0, index.localDissimilarity(seedNode, cluster, groundTruth))
			<< "The Jaccard index has a value of 0";

	cluster = {0,1,2,3};
	EXPECT_EQ(0.4, index.localDissimilarity(seedNode, cluster, groundTruth))
			<< "The Jaccard index has a value of 0.4";
}

TEST_F(SCDGTest, JaccardTest) {
	std::unordered_map<node, std::unordered_set<node>> communities;
	communities.insert( { 1, { 1, 4 } });
	communities.insert( { 2, { 2, 5 } });
	Clustering groundTruth(8);
	JaccardIndex index;

	groundTruth.addToCluster(0, 1);
	groundTruth.addToCluster(0, 6);
	groundTruth.addToCluster(1, 2);
	groundTruth.addToCluster(1, 4);
	groundTruth.addToCluster(1, 7);

	EXPECT_LE(0.33333, index.getDissimilarity(communities, groundTruth))
			<< "The Jaccard index has a value of 1/3";
	EXPECT_GE(0.33334, index.getDissimilarity(communities, groundTruth))
			<< "The Jaccard index has a value of 1/3";
}

TEST_F(SCDGTest, localPrecisionTest) {
	node seedNode = 0;
	std::unordered_set<node> cluster = { 0, 1, 2 };
	Clustering groundTruth(5);
	Precision index;

	groundTruth.addToCluster(0, 0);
	groundTruth.addToCluster(0, 1);
	groundTruth.addToCluster(0, 3);
	groundTruth.addToCluster(1, 4);

	EXPECT_LE(0.66666, index.localDissimilarity(seedNode, cluster, groundTruth))
			<< "The Precision has a value of 2/3";
	EXPECT_GE(0.66667, index.localDissimilarity(seedNode, cluster, groundTruth))
			<< "The Precision has a value of 2/3";
}

TEST_F(SCDGTest, PrecisionTest) {
	std::unordered_map<node, std::unordered_set<node>> communities;
	communities.insert( { 1, { 1, 4 } });
	communities.insert( { 2, { 2, 5 } });
	Clustering groundTruth(8);
	Precision index;

	groundTruth.addToCluster(0, 1);
	groundTruth.addToCluster(0, 6);
	groundTruth.addToCluster(1, 2);
	groundTruth.addToCluster(1, 4);
	groundTruth.addToCluster(1, 7);

	EXPECT_EQ(0.5, index.getDissimilarity(communities, groundTruth))
			<< "The Precision has a value of 0.5";
}

TEST_F(SCDGTest, localRecallTest) {
	node seedNode = 0;
	std::unordered_set<node> cluster = { 0, 1, 2 };
	Clustering groundTruth(5);
	Recall index;

	groundTruth.addToCluster(0, 0);
	groundTruth.addToCluster(0, 1);
	groundTruth.addToCluster(0, 3);
	groundTruth.addToCluster(1, 4);

	EXPECT_LE(0.66666, index.localDissimilarity(seedNode, cluster, groundTruth))
			<< "The Recall has a value of 2/3";
	EXPECT_GE(0.66667, index.localDissimilarity(seedNode, cluster, groundTruth))
			<< "The Recall has a value of 2/3";
}

TEST_F(SCDGTest, RecallTest) {
	std::unordered_map<node, std::unordered_set<node>> communities;
	communities.insert( { 1, { 1, 4 } });
	communities.insert( { 2, { 2, 5 } });
	Clustering groundTruth(8);
	Recall index;

	groundTruth.addToCluster(0, 1);
	groundTruth.addToCluster(0, 6);
	groundTruth.addToCluster(1, 2);
	groundTruth.addToCluster(1, 4);
	groundTruth.addToCluster(1, 7);

	EXPECT_EQ(0.4, index.getDissimilarity(communities, groundTruth))
			<< "The Precision has a value of 0.4";
}

// Test Acceptability
TEST_F(SCDGTest, testNodeClusterSimilarity) {

	Graph G(5);
	G.addEdge(0, 0);
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(1, 2);
	G.addEdge(1, 4);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 4);

	std::unordered_set<node> first = { 0 };
	std::unordered_set<node> shell_1 = { 1, 2, 3 };
	NodeClusterSimilarity ncs1(G, first, shell_1);
	double ncsOne = ncs1.getValue(1);
	EXPECT_EQ(0.6, ncsOne) << "1-clustering should have similarity of 0.6";

	std::unordered_set<node> second = { 1 };
	std::unordered_set<node> shell_2 = { 2, 0, 4 };
	NodeClusterSimilarity ncs2(G, second, shell_2);
	double ncsTwo = ncs2.getValue(0);
	EXPECT_EQ(0.6, ncsTwo) << "2-clustering should have similarity of 0.6";

	std::unordered_set<node> third = { 0, 1, 2 };
	std::unordered_set<node> shell_3 = { 3, 4 };
	NodeClusterSimilarity ncs3(G, third, shell_3);
	double ncsThree = ncs3.getValue(3);
	EXPECT_GE(0.8, ncsThree) << "3-clustering should have similarity of 0.8";

	std::unordered_set<node> fourth = { 0, 1, 2, 3 };
	std::unordered_set<node> shell_4 = { 4 };
	NodeClusterSimilarity ncs4(G, fourth, shell_4);
	double ncsFour = ncs4.getValue(4);
	EXPECT_EQ(0.8, ncsFour) << "4-clustering should have similarity of 0.8";

	std::unordered_set<node> fifth = { 0, 1 };
	std::unordered_set<node> shell_5 = { 2, 3, 4 };
	NodeClusterSimilarity ncs5(G, fifth, shell_5);
	double ncsFive = ncs5.getValue(2);
	EXPECT_EQ(1, ncsFive) << "5-clustering should have similarity of 1";
}
// Test Qualilty Measures
TEST_F(SCDGTest, testLocalModularity) {
	Graph G(5);
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(0, 4);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(1, 4);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 4);
	LocalModularity mod(G);
	std::unordered_set<node> set;
	set.insert(0);
	EXPECT_EQ(0, mod.getQuality(set))
				<< "The community should have a local modularity of 0";
	set.erase(0);
	set.insert(4);
	EXPECT_EQ(0.25, mod.getQuality(set))
				<< "The community should have a local modularity of 0.25";

	set.erase(4);
	set.insert(0);
	set.insert(1);
	EXPECT_GE(0.16667, mod.getQuality(set))
			<< "The community should have a local modularity of 1/6";
	EXPECT_LE(0.16666, mod.getQuality(set))
			<< "The community should have a local modularity of 1/6";
	set.erase(1);
	set.insert(4);
	EXPECT_GE(0.33334, mod.getQuality(set))
			<< "The community should have a local modularity of 1/3";
	EXPECT_LE(0.33333, mod.getQuality(set))
			<< "The community should have a local modularity of 1/3";

	set.erase(4);
	set.insert(1);
	set.insert(2);
	EXPECT_EQ(0.5, mod.getQuality(set))
			<< "The community should have a local modularity of 0.5";
	set.erase(2);
	set.insert(4);
	EXPECT_GE(0.66667, mod.getQuality(set))
			<< "The community should have a local modularity of 2/3";
	EXPECT_LE(0.66666, mod.getQuality(set))
			<< "The community should have a local modularity of 2/3";

	set.erase(4);
	set.insert(2);
	set.insert(3);
	EXPECT_EQ(1.5, mod.getQuality(set))
			<< "The community should have a local modularity of 1.5";
	set.erase(3);
	set.insert(4);
	EXPECT_EQ(1.75, mod.getQuality(set))
			<< "The community should have a local modularity of 1.75";


	Graph G1(6);
		G1.addEdge(0, 1);
		G1.addEdge(0, 2);
		G1.addEdge(1, 2);
		G1.addEdge(1, 3);
		G1.addEdge(2, 3);
		G1.addEdge(2, 4);
		G1.addEdge(2, 5);
		G1.addEdge(4, 5);
		std::unordered_set<node> set1 = {0, 1, 2, 3};
		LocalModularity mod1(G1);
		EXPECT_EQ(0.625, mod1.getQuality(set1))
				<< "The community should have a local modularity of 5/8";
		set1.erase(3);
		set1.insert(4);
		EXPECT_EQ(0.75, mod1.getQuality(set1))
				<< "The community should have a local modularity of 0.75";
}

// Test QUalityObjective
TEST_F(SCDGTest, testLocalMdularityM) {

	Graph G(5);
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(0, 4);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(1, 4);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 4);

	std::unordered_set<node> community = { };
	std::unordered_map<node, count> bound = {};
	LocalModularityM mod(G, community, bound);
	mod.degSum = 11;
	mod.volume = 0;
	mod.nInternEdges = 0;
	EXPECT_EQ(0, mod.getValue(0)[0])
			<< "The community should have a local modularity of 0";
	EXPECT_EQ(0.25, mod.getValue(4)[0])
			<< "The community should have a local modularity of 0.25";

	community.insert(0);
	mod.nBoundaryEdges = 4;
	mod.volume = 4;
	mod.nInternEdges = 0;

	EXPECT_EQ(0, mod.getValue(0)[0])
			<< "The community should have a local modularity of 0";
	EXPECT_GE(0.16667, mod.getValue(1)[0])
			<< "The community should have a local modularity of 1/6";
	EXPECT_LE(0.16666, mod.getValue(1)[0])
			<< "The community should have a local modularity of 1/6";
	EXPECT_GE(0.33334, mod.getValue(4)[0])
			<< "The community should have a local modularity of 1/3";
	EXPECT_LE(0.33333, mod.getValue(4)[0])
			<< "The community should have a local modularity of 1/3";

	community.insert(1);
	mod.nBoundaryEdges = 6;
	mod.volume = 8;
	mod.nInternEdges = 1;
	EXPECT_EQ(0.5, mod.getValue(2)[0])
			<< "The community should have a local modularity of 0.5";
	EXPECT_GE(0.66667, mod.getValue(4)[0])
			<< "The community should have a local modularity of 2/3";
	EXPECT_LE(0.66666, mod.getValue(4)[0])
			<< "The community should have a local modularity of 2/3";

	community.insert(2);
	mod.nBoundaryEdges = 6;
	mod.volume = 12;
	mod.nInternEdges = 3;
	EXPECT_EQ(1.5, mod.getValue(3)[0])
			<< "The community should have a local modularity of 1.5";
	EXPECT_EQ(1.75, mod.getValue(4)[0])
			<< "The community should have a local modularity of 1.75";

	community.insert(3);
	mod.nBoundaryEdges = 4;
	mod.volume = 16;
	mod.nInternEdges = 6;
	EXPECT_EQ(11, mod.getValue(4)[0])
			<< "The community should have a local modularity of 11";
}

TEST_F(SCDGTest, testLocalMdularityL) {

	Graph G(5);
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(0, 4);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(1, 4);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 4);

	std::unordered_set<node> community = { };
	std::unordered_map<node, count> bound = {};
	LocalModularityL mod(G, community, bound);
	mod.degSum = 11;
	mod.volume = 0;
	mod.nInternEdges = 0;
	EXPECT_EQ(0, mod.getValue(0)[0])
			<< "The community should have a local modularity of 0";

	EXPECT_EQ(0.25, mod.getValue(4)[0])
			<< "The community should have a local modularity of 0.25";

	community.insert(0);
	bound = {{0, 4}};
	mod.nBoundaryEdges = 4;
	mod.volume = 4;
	mod.nInternEdges = 0;

	EXPECT_EQ(0, mod.getValue(0)[0])
			<< "The community should have a local modularity of 0";
	EXPECT_GE(0.16667, mod.getValue(1)[0])
			<< "The community should have a local modularity of 1/6";
	EXPECT_LE(0.16666, mod.getValue(1)[0])
			<< "The community should have a local modularity of 1/6";
	EXPECT_GE(0.33334, mod.getValue(4)[0])
			<< "The community should have a local modularity of 1/3";
	EXPECT_LE(0.33333, mod.getValue(4)[0])
			<< "The community should have a local modularity of 1/3";

	community.insert(1);
	bound = {{0, 3}, {1, 3}};
	mod.nBoundaryEdges = 6;
	mod.volume = 8;
	mod.nInternEdges = 1;
	EXPECT_EQ(0.5, mod.getValue(2)[0])
			<< "The community should have a local modularity of 0.5";
	EXPECT_GE(0.66667, mod.getValue(4)[0])
			<< "The community should have a local modularity of 2/3";
	EXPECT_LE(0.66666, mod.getValue(4)[0])
			<< "The community should have a local modularity of 2/3";

	community.insert(2);
	bound = {{0, 2}, {1, 2}, {2, 2}};
	mod.nBoundaryEdges = 6;
	mod.volume = 12;
	mod.nInternEdges = 3;
	EXPECT_EQ(1.5, mod.getValue(3)[0])
			<< "The community should have a local modularity of 1.5";
	EXPECT_EQ(1.75, mod.getValue(4)[0])
			<< "The community should have a local modularity of 1.75";

	community.insert(3);
	bound = {{0, 1}, {1, 1},{2, 1}, {3, 1}};
	mod.nBoundaryEdges = 4;
	mod.volume = 16;
	mod.nInternEdges = 6;
	EXPECT_EQ(11, mod.getValue(4)[0])
			<< "The community should have a local modularity of 11";

	Graph G1(6);
	G1.addEdge(0, 1);
	G1.addEdge(0, 2);
	G1.addEdge(1, 2);
	G1.addEdge(1, 3);
	G1.addEdge(2, 3);
	G1.addEdge(2, 4);
	G1.addEdge(2, 5);
	G1.addEdge(4, 5);
	std::unordered_set<node> community1 = {0, 1, 2};
	std::unordered_map<node,count> bound1 = {{1, 1},{2, 3}};
	LocalModularityL mod1(G1, community1, bound1);
	mod1.nBoundaryEdges = 4;
	mod1.volume = 10;
	mod1.nInternEdges = 3;
	EXPECT_EQ(0.625, mod1.getValue(3)[0])
			<< "The community should have a local modularity of 5/8";
	EXPECT_EQ(0.75, mod1.getValue(4)[0])
				<< "The community should have a local modularity of 0.75";
}

TEST_F(SCDGTest, testConductance) {

	Graph G(5);
	G.addEdge(0, 0);
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(1, 2);
	G.addEdge(1, 4);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 4);

	std::unordered_set<node> first, second, third, fourth, fifth;
	first = {};
	std::unordered_map<node, count> bound1 = {};
	Conductance conductance1(G, first, bound1);
	second = {0};
	std::unordered_map<node, count> bound2 = {{0, 3}};
	Conductance conductance2(G, second, bound2);
	conductance2.volume = 4;
	conductance2.nBoundaryEdges = 3;
	third = {0,1};
	std::unordered_map<node, count> bound3 = {{0, 2},{1, 2}};
	Conductance conductance3(G, third, bound3);
	conductance3.volume = 7;
	conductance3.nBoundaryEdges = 4;
	fourth = {0,1,2};
	std::unordered_map<node, count> bound4 = {{0, 1},{1, 1},{2, 2}};
	Conductance conductance4(G, fourth, bound4);
	conductance4.volume = 11;
	conductance4.nBoundaryEdges = 4;
	fifth = {0,1,2,3};
	std::unordered_map<node, count> bound5 = {{1, 1},{2, 1},{3, 1}};
	Conductance conductance5(G, fifth, bound5);
	conductance5.volume = 14;
	conductance5.nBoundaryEdges = 3;

	double condOne = conductance1.getValue(0)[0];
	double condTwo = conductance2.getValue(1)[0];
	double condThree = conductance3.getValue(2)[0];
	double condFour = conductance4.getValue(3)[0];
	double condFive = conductance5.getValue(4)[0];

	EXPECT_EQ(0.75, 1-condOne)
			<< "1-clustering should have conductance of 0.75";

	EXPECT_GE(0.571429, 1-condTwo)
			<< "2-clustering should have conductance of 4/7";
	EXPECT_LE(0.571428, 1-condTwo)
			<< "2-clustering should have conductance of 4/7";

	EXPECT_GE(0.666667, 1-condThree)
			<< "3-clustering should have conductance of 2/3";
	EXPECT_LE(0.666666, 1-condThree)
			<< "3-clustering should have conductance of 2/3";
	EXPECT_EQ(1, 1-condFour) << "4-clustering should have conductance of 1";

	EXPECT_EQ(1, 1-condFive) << "5-clustering should have conductance of 1";

}

// Test Trimming

TEST_F(SCDGTest, testBoundarySharpnessTrimming) {

	Graph G(11);
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(0, 4);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(1, 5);
	G.addEdge(2, 4);
	G.addEdge(3, 10);
	G.addEdge(4, 6);
	G.addEdge(4, 7);
	G.addEdge(5, 8);
	G.addEdge(5, 9);

	std::unordered_set<node> cluster = { 0, 1, 2, 3, 4, 5 };
	DummyTrimming trim;
	std::unordered_set<node> community = trim.run(cluster, G);
	EXPECT_EQ(6, community.size()) << "The community has 6 nodes";

	BoundarySharpness trim2;
	std::unordered_set<node> community2 = trim2.run(cluster, G);
	EXPECT_EQ(5, community2.size()) << "The community has 5 nodes";

}

// Test Algorithms

//Test with dummy acceptance and conductance and dummy trimming
TEST_F(SCDGTest, testGreedyCommunityExpansion) {

	Graph G(12);
	// add clique
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(2, 3);
	std::unordered_set<node> tmp = {};
	std::unordered_map<node, count> bound =  {};
	DummySimilarity similarity(G, tmp, tmp);
	Conductance objective(G, tmp, bound);
	DummyTrimming trimming;
	GreedyCommunityExpansion GCE(G, similarity, objective, trimming);

	std::unordered_set<node> community = GCE.expandSeed(0);
	EXPECT_EQ(2, community.size()) << "The community should have 2 nodes";

	//add satelites
	G.addEdge(0, 4);
	G.addEdge(1, 5);
	G.addEdge(2, 6);
	G.addEdge(3, 7);
	DummySimilarity similarity1(G, tmp, tmp);
	Conductance objective1(G, tmp, bound);
	GreedyCommunityExpansion GCE1(G, similarity1, objective1, trimming);

	community = GCE1.expandSeed(0);
	EXPECT_EQ(4, community.size()) << "The community should have 4 nodes";
	community = GCE1.expandSeed(6);
	EXPECT_EQ(4, community.size()) << "The community should have 4 nodes";

	// add another clique
	G.addEdge(6, 8);
	G.addEdge(8, 9);
	G.addEdge(8, 10);
	G.addEdge(8, 11);
	G.addEdge(9, 10);
	G.addEdge(9, 11);
	G.addEdge(10, 11);
	DummySimilarity similarity2(G, tmp, tmp);
	Conductance objective2(G, tmp, bound);
	GreedyCommunityExpansion GCE2(G, similarity2, objective2, trimming);
	community = GCE2.expandSeed(0);
	EXPECT_EQ(7, community.size()) << "The community should have 7 nodes";
	community = GCE2.expandSeed(4);
	EXPECT_EQ(7, community.size()) << "The community should have 7 nodes";
	community = GCE2.expandSeed(6);
	EXPECT_EQ(8, community.size()) << "The community should have 8 nodes";
	community = GCE2.expandSeed(8);
	EXPECT_EQ(5, community.size()) << "The community should have 5 nodes";
	community = GCE2.expandSeed(9);
	EXPECT_EQ(5, community.size()) << "The community should have 5 nodes";
}

TEST_F(SCDGTest, testTSelectiveSCAN) {
	Graph G(10);
	// add clique
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(0, 4);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(1, 4);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 5);
	G.addEdge(4, 6);
	G.addEdge(4, 7);
	G.addEdge(4, 8);
	G.addEdge(5, 6);
	G.addEdge(5, 7);
	G.addEdge(5, 8);
	G.addEdge(6, 7);
	G.addEdge(6, 8);
	G.addEdge(7, 8);
	G.addEdge(0, 9);

	std::unordered_set<node> seeds = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	NeighborhoodDistance distMeasure(G);
	SelectiveSCAN GCE(G, distMeasure);
	std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> result = GCE.run(seeds);
	std::unordered_set<node> set1 = { 0, 1, 2, 3 };
	std::unordered_set<node> set2 = { 5, 6, 7, 8 };

	EXPECT_EQ(0, result.find(4)->second.first.size())
			<< "The node 4 has no community";
	EXPECT_EQ(0, result.find(9)->second.first.size())
			<< "The node 9 has no community";

	for (int i = 0; i < 4; i++) {
		EXPECT_EQ(4, result.find(i)->second.first.size()) << "The community of node"
				<< i << " has 4 nodes";
		for (int j = 0; j < 4; j++) {
			EXPECT_TRUE(result.find(i)->second.first.find(j) != result.find(i)->second.first.end())
					<< "The node" << j << " belongs to the community of node"
					<< i;
		}
	}
	for (int i = 5; i < 9; i++) {
		EXPECT_EQ(4, result.find(i)->second.first.size()) << "The community of node"
				<< i << "has 4 nodes";
		for (int j = 5; j < 9; j++) {
			EXPECT_TRUE(result.find(i)->second.first.find(j) != result.find(i)->second.first.end())
					<< "The node" << j << "belongs to the community of node"
					<< i;
		}
	}
}


TEST_F(SCDGTest, testTSelectiveSCANwithTemplates) {
	Graph G(10);
	// add clique
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(0, 4);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(1, 4);
	G.addEdge(2, 3);
	G.addEdge(2, 4);
	G.addEdge(3, 4);
	G.addEdge(4, 5);
	G.addEdge(4, 6);
	G.addEdge(4, 7);
	G.addEdge(4, 8);
	G.addEdge(5, 6);
	G.addEdge(5, 7);
	G.addEdge(5, 8);
	G.addEdge(6, 7);
	G.addEdge(6, 8);
	G.addEdge(7, 8);
	G.addEdge(0, 9);

	std::unordered_set<node> seeds = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

	Parameters param;
	TSelectiveSCAN<TNeighborhoodDistance> GCE(G, param);
	std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> result = GCE.run(seeds);
	std::unordered_set<node> set1 = { 0, 1, 2, 3 };
	std::unordered_set<node> set2 = { 5, 6, 7, 8 };

	EXPECT_EQ(0, result.find(4)->second.first.size())
			<< "The node 4 has no community";
	EXPECT_EQ(0, result.find(9)->second.first.size())
			<< "The node 9 has no community";

	for (int i = 0; i < 4; i++) {
		EXPECT_EQ(4, result.find(i)->second.first.size()) << "The community of node"
				<< i << "has 4 nodes";
		for (int j = 0; j < 4; j++) {
			EXPECT_TRUE(result.find(i)->second.first.find(j) != result.find(i)->second.first.end())
					<< "The node" << j << "belongs to the community of node"
					<< i;
		}
	}
	for (int i = 5; i < 9; i++) {
		EXPECT_EQ(4, result.find(i)->second.first.size()) << "The community of node"
				<< i << "has 4 nodes";
		for (int j = 5; j < 9; j++) {
			EXPECT_TRUE(result.find(i)->second.first.find(j) != result.find(i)->second.first.end())
					<< "The node" << j << "belongs to the community of node"
					<< i;
		}
	}
}

TEST_F(SCDGTest, testTGreedyCommunityExpansionWithTemplates) {

	Graph G(12);

	// add clique
	G.addEdge(0, 1);
	G.addEdge(0, 2);
	G.addEdge(0, 3);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(2, 3);

	TGreedyCommunityExpansion<TConductance, TDummyAcceptability, DummyTrimming> GCE1(
			G);
	std::unordered_set<node> community = GCE1.expandSeed(0);
	EXPECT_EQ(2, community.size()) << "The community should have 2 nodes";

	G.addEdge(0, 4);
	G.addEdge(1, 5);
	G.addEdge(2, 6);
	G.addEdge(3, 7);

	TGreedyCommunityExpansion<TConductance, TDummyAcceptability, DummyTrimming> GCE2(
			G);
	community = GCE2.expandSeed(0);
	EXPECT_EQ(4, community.size()) << "The community should have 4 nodes";
	community = GCE2.expandSeed(6);
	EXPECT_EQ(4, community.size()) << "The community should have 4 nodes";

//	add another clique
	G.addEdge(6, 8);
	G.addEdge(8, 9);
	G.addEdge(8, 10);
	G.addEdge(8, 11);
	G.addEdge(9, 10);
	G.addEdge(9, 11);
	G.addEdge(10, 11);

	TGreedyCommunityExpansion<TConductance, TDummyAcceptability, DummyTrimming> GCE3(
			G);
	community = GCE3.expandSeed(0);
	EXPECT_EQ(7, community.size()) << "The community should have 7 nodes";
	community = GCE3.expandSeed(4);
	EXPECT_EQ(7, community.size()) << "The community should have 7 nodes";
	community = GCE3.expandSeed(6);
	EXPECT_EQ(8, community.size()) << "The community should have 8 nodes";
	community = GCE3.expandSeed(8);
	EXPECT_EQ(5, community.size()) << "The community should have 5 nodes";
	community = GCE3.expandSeed(9);
	EXPECT_EQ(5, community.size()) << "The community should have 5 nodes";
}

TEST_F(SCDGTest, lol) {
	Graph G;
	std::string path = "input/graph1000.dat";
	std::string path1 = "input/graph1000.graph";
	EdgeListReader graphReader;
	G = graphReader.read(path);
	METISGraphWriter writer;
	writer.write(G, path1);
}

//TEST_F(SCDGTest, trySelectiveSCANWithSeedSet) {
//
//	METISGraphReader reader;
//	Graph G = reader.read("input/caidaRouterLevel.graph");
//	Parameters param;
//	TSelectiveSCAN<TNeighborhoodDistance> GCE(G, param);
//	RandomSeedSet randSeeds(G);
//	for (int i = 0; i < 100; i++) {
//		std::unordered_set<node> seeds = randSeeds.getSeeds(77);
//		std::unordered_map<node, std::unordered_set<node>> result = GCE.run(
//				seeds);
//	}
//}
//
//TEST_F(SCDGTest, trySelectiveSCANWithSeedSetAndTemplates) {
//
//	METISGraphReader reader;
//	Graph G = reader.read("input/pgp.graph");
//	AlgebraicDistance distMeasure(G, 10, 30, 0.5, 2);
//	SelectiveSCAN GCE(G, distMeasure, 0.005, 3);
//
//	RandomSeedSet randSeeds(G);
//
//		std::unordered_set<node> seeds = {7327};//randSeeds.getSeeds(100);
//		std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> result = GCE.run(
//				seeds);
//		for(auto u:result) {
//			std::cout<<u.second.first.size()<<std::endl;
//		}
//}

//}
//TEST_F(SCDGTest, trySelectlates) {
//
//	METISGraphReader reader;
//	Graph G = reader.read("input/pgp.graph");
//	Parameters param;
//	Conduct x(G);
//	TGreedyCommunityExpansion<TConductance, DummySimilarity, DummyTrimming> GCE(G);
//	RandomSeedSet randSeeds(G);
//	//for(int i = 0; i < 1000; i++) {
//		std::unordered_set<node> seeds = {7592,5196,6239,4133,9648,6975,3822};//randSeeds.getSeeds(100);
//		std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> result = GCE.run(
//				seeds);
//for(auto u: result){
//	std::cout<< x.getQuality(u.second.first)<<std::endl;
//	}
//}


} /* namespace NetworKit */

#endif /*NOGTEST */
