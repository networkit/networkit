/*
 * MatcherGTest.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#include "MatcherGTest.h"
#include "../Matcher.h"
#include "../Matching.h"
#include "../PathGrowingMatcher.h"
#include "../LocalMaxMatcher.h"
#include "../../graph/Graph.h"
#include "../../io/DibapGraphReader.h"
#include "../../io/METISGraphReader.h"


namespace NetworKit {


TEST_F(MatcherGTest, testLocalMaxMatching) {
	count n = 50;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u,v);
	});

	LocalMaxMatcher localMaxMatcher(G);

	TRACE("Start localMax matching");
	Matching M = localMaxMatcher.run();
	TRACE("Finished localMax matching");

	count numExpEdges = n / 2;
	bool isProper = M.isProper(G);
	EXPECT_TRUE(isProper);
	EXPECT_EQ(M.size(), numExpEdges);

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
	DibapGraphReader reader;
	Graph airfoil1 = reader.read("input/airfoil1.gi");
	LocalMaxMatcher lmm(airfoil1);
	M = lmm.run();
	isProper = M.isProper(airfoil1);
	EXPECT_TRUE(isProper);
	DEBUG("LocalMax on airfoil1 produces matching of size: " , M.size());
#endif
}


TEST_F(MatcherGTest, testPgaMatching) {
	count n = 50;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u,v);
	});
	PathGrowingMatcher pgaMatcher(G);

	DEBUG("Start PGA matching on 50-clique");
	Matching M = pgaMatcher.run();

	count numExpEdges = n / 2;
	bool isProper = M.isProper(G);
	EXPECT_TRUE(isProper);
	EXPECT_EQ(M.size(), numExpEdges);
	DEBUG("Finished PGA matching on 50-clique");


#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
	DibapGraphReader reader;
	Graph airfoil1 = reader.read("input/airfoil1.gi");
	PathGrowingMatcher pga2(airfoil1);
	M = pga2.run();
	isProper = M.isProper(airfoil1);
	EXPECT_TRUE(isProper);
	DEBUG("PGA on airfoil1 produces matching of size: " , M.size());
#endif
}

TEST_F(MatcherGTest, tryValidMatching) {
	METISGraphReader reader;
	Graph G = reader.read("coAuthorsDBLP.graph");

	LocalMaxMatcher pmatcher(G);
	Matching M = pmatcher.run();

	bool isProper = M.isProper(G);
	EXPECT_TRUE(isProper);
}


} // namespace EnsembleClustering

#endif
