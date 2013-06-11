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
//	Conductance conductance;
//	std::unordered_set<node> first, second, third, fourth;
//	first = {0};
//	second = {0,1};
//	third = {0,1,2,3,4};
//	fourth = {0,1,2,3};
//
//	double condThree = conductance.getQuality(third, G);
//
//	EXPECT_EQ(1.0, condThree) << "1-clustering should have modularity of 1.0";
//
//}

} /* namespace NetworKit */

#endif /* NOGTEST */
