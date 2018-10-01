/*
 * MatcherGTest.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Henning
 */

#include <gtest/gtest.h>

#include "../Matcher.h"
#include "../Matching.h"
#include "../PathGrowingMatcher.h"
#include "../LocalMaxMatcher.h"
#include "../../graph/Graph.h"
#include "../../io/DibapGraphReader.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Random.h"

namespace NetworKit {

class MatcherGTest: public testing::Test {};

TEST_F(MatcherGTest, testLocalMaxMatching) {
	count n = 50;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u,v);
	});

	LocalMaxMatcher localMaxMatcher(G);

	TRACE("Start localMax matching");
	localMaxMatcher.run();
	Matching M = localMaxMatcher.getMatching();
	TRACE("Finished localMax matching");

	count numExpEdges = n / 2;
	bool isProper = M.isProper(G);
	EXPECT_TRUE(isProper);
	EXPECT_EQ(M.size(G), numExpEdges);

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
	DibapGraphReader reader;
	Graph airfoil1 = reader.read("input/airfoil1.gi");
	LocalMaxMatcher lmm(airfoil1);
	lmm.run();
	M = lmm.getMatching();
	isProper = M.isProper(airfoil1);
	EXPECT_TRUE(isProper);
	DEBUG("LocalMax on airfoil1 produces matching of size: " , M.size(G));
#endif
}

TEST_F(MatcherGTest, testLocalMaxMatchingDirectedWarning) {
	Graph G(2, false, true);
	G.addEdge(0,1);
	EXPECT_THROW(LocalMaxMatcher localMaxMatcher(G), std::runtime_error);
}

TEST_F(MatcherGTest, testPgaMatchingOnWeightedGraph) {
	count n = 50;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u,v, Aux::Random::real());
	});
	PathGrowingMatcher pgaMatcher(G);
	EXPECT_NO_THROW(pgaMatcher.run());
}

TEST_F(MatcherGTest, testPgaMatchingWithSelfLoops) {
	count n = 50;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u,v, Aux::Random::real());
	});
	G.forNodes([&](node u){
		G.addEdge(u,u);
	});
	EXPECT_THROW(PathGrowingMatcher pgaMatcher(G),std::invalid_argument);
}


TEST_F(MatcherGTest, testPgaMatching) {
	count n = 50;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u,v);
	});
	PathGrowingMatcher pgaMatcher(G);

	DEBUG("Start PGA matching on 50-clique");

	pgaMatcher.run();
	Matching M = pgaMatcher.getMatching();

	count numExpEdges = n / 2;
	bool isProper = M.isProper(G);
	EXPECT_TRUE(isProper);
	EXPECT_EQ(M.size(G), numExpEdges);
	DEBUG("Finished PGA matching on 50-clique");


#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
	DibapGraphReader reader;
	Graph airfoil1 = reader.read("input/airfoil1.gi");
	PathGrowingMatcher pga2(airfoil1);
	pga2.run();
	M = pga2.getMatching();
	isProper = M.isProper(airfoil1);
	EXPECT_TRUE(isProper);
	DEBUG("PGA on airfoil1 produces matching of size: " , M.size(G));
#endif
}

TEST_F(MatcherGTest, debugValidMatching) {
	METISGraphReader reader;
	Graph G = reader.read("coAuthorsDBLP.graph");

	LocalMaxMatcher pmatcher(G);
	pmatcher.run();
	Matching M = pmatcher.getMatching();

	bool isProper = M.isProper(G);
	EXPECT_TRUE(isProper);
}

} // namespace NetworKit
