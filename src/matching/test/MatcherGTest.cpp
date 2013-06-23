/*
 * MatcherGTest.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#include "MatcherGTest.h"


namespace NetworKit {


MatcherGTest::MatcherGTest() {

}

MatcherGTest::~MatcherGTest() {

}


TEST_F(MatcherGTest, testLocalMaxMatching) {
	count n = 50;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);
	LocalMaxMatcher localMaxMatcher(none);

	TRACE("Start localMax matching");
	Matching M = localMaxMatcher.run(G);
	TRACE("Finished localMax matching");

	count numExpEdges = n / 2;
	bool isProper = M.isProper(G);
	EXPECT_TRUE(isProper);
	EXPECT_EQ(M.size(), numExpEdges);

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
	DibapGraphReader reader;
	Graph airfoil1 = reader.read("input/airfoil1.gi");
	M = localMaxMatcher.run(airfoil1);
	isProper = M.isProper(airfoil1);
	EXPECT_TRUE(isProper);
	DEBUG("LocalMax on airfoil1 produces matching of size: " << M.size());
#endif
}


TEST_F(MatcherGTest, testPgaMatching) {
	count n = 50;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);
	PathGrowingMatcher pgaMatcher;

	TRACE("Start PGA matching");
	Matching M = pgaMatcher.run(G);
	TRACE("Finished PGA matching");

	count numExpEdges = n / 2;
	bool isProper = M.isProper(G);
	EXPECT_TRUE(isProper);
	EXPECT_EQ(M.size(), numExpEdges);

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
	DibapGraphReader reader;
	Graph airfoil1 = reader.read("input/airfoil1.gi");
	M = pgaMatcher.run(airfoil1);
	isProper = M.isProper(airfoil1);
	EXPECT_TRUE(isProper);
	DEBUG("PGA on airfoil1 produces matching of size: " << M.size());
#endif
}



} // namespace EnsembleClustering

#endif
