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
	GreedyCommunityExpansion::Conductance conductance1(G, &first);
	second = {0};
	GreedyCommunityExpansion::Conductance conductance2(G, &second);
	third = {0,1};
	GreedyCommunityExpansion::Conductance conductance3(G, &third);
	fourth = {0,1,2};
	GreedyCommunityExpansion::Conductance conductance4(G, &fourth);
	fifth = {0,1,2,3};
	GreedyCommunityExpansion::Conductance conductance5(G, &fifth);

	double condOne = conductance1.getValue(0);
	double condTwo = conductance2.getValue(1);
	double condThree = conductance3.getValue(2);
	double condFour = conductance4.getValue(3);
	double condFive = conductance5.getValue(4);


	EXPECT_EQ(0.75, condOne) << "1-clustering should have conductance of 0.75";

	EXPECT_GE(0.571429, condTwo) << "2-clustering should have conductance of 4/7";
	EXPECT_LE(0.571428, condTwo) << "2-clustering should have conductance of 4/7";

	EXPECT_GE(0.666667, condThree) << "3-clustering should have conductance of 2/3";
	EXPECT_LE(0.666666, condThree) << "3-clustering should have conductance of 2/3";

	EXPECT_EQ(1, condFour) << "4-clustering should have conductance of 1";

	EXPECT_EQ(1, condFive) << "5-clustering should have conductance of 1";

}

SCDGTest::SCDGTest() {
}

SCDGTest::~SCDGTest() {
}
} /* namespace NetworKit */

#endif /*NOGTEST */
