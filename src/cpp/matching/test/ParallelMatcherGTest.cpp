/*
 * ParallelMatcherGTest.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#include "ParallelMatcherGTest.h"

#include "../ParallelMatcher.h"
#include "../Matching.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {


ParallelMatcherGTest::ParallelMatcherGTest() {
	// TODO Auto-generated constructor stub

}

ParallelMatcherGTest::~ParallelMatcherGTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(ParallelMatcherGTest, tryValidMatching) {
	METISGraphReader reader;
	Graph G = reader.read("coAuthorsDBLP.graph");

	LocalMaxMatcher pmatcher(none);
	Matching M = pmatcher.run(G);

	bool isProper = M.isProper(G);
	EXPECT_TRUE(isProper);
}



} // namespace EnsembleClustering


#endif
