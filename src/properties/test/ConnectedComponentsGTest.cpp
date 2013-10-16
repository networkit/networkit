/*
 * ConnectedComponentsGTest.cpp
 *
 *  Created on: Sep 16, 2013
 *      Author: Maximilian Vogel
 */
#ifndef NOGTEST

#include "ConnectedComponentsGTest.h"

namespace NetworKit {

ConnectedComponentsGTest::ConnectedComponentsGTest() {
	// TODO Auto-generated constructor stub
}

ConnectedComponentsGTest::~ConnectedComponentsGTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(ConnectedComponentsGTest, TestRunByCall) {
	// construct graph
	Graph g;
	for (count i = 0; i < 20; i++) {
		g.addNode();
	}
	g.addEdge(0,1,0);
	g.addEdge(1,2,0);
	g.addEdge(2,4,0);
	g.addEdge(4,8,0);
	g.addEdge(8,16,0);
	g.addEdge(16,19,0);

	g.addEdge(3,5,0);
	g.addEdge(5,6,0);
	g.addEdge(6,7,0);
	g.addEdge(7,9,0);

	g.addEdge(10,11,0);
	g.addEdge(10,18,0);
	g.addEdge(10,12,0);
	g.addEdge(18,17,0);

	g.addEdge(13,14,0);

	// initialize ConnectedComponents
	ConnectedComponents ccs;
	ccs.run(g);

	// check result
	EXPECT_TRUE(ccs.numberOfComponents() == 5);
	EXPECT_TRUE(ccs.componentOfNode(0) == ccs.componentOfNode(19));
	EXPECT_TRUE(ccs.componentOfNode(3) == ccs.componentOfNode(7));
}


TEST_F(ConnectedComponentsGTest, TestGraphWithoutEdges) {
	// construct graph
	Graph g;
	for (count i = 0; i < 20; i++) {
		g.addNode();
	}

	// initialize ConnectedComponents
	ConnectedComponents ccs;
	ccs.run(g);

	// check result
	EXPECT_TRUE(ccs.numberOfComponents() == 20);
	EXPECT_TRUE(ccs.componentOfNode(0) != ccs.componentOfNode(19));
	EXPECT_TRUE(ccs.componentOfNode(3) != ccs.componentOfNode(7));
}

TEST_F(ConnectedComponentsGTest, TestGetComponent) {
	// construct graph
	Graph g;
	for (count i = 0; i < 20; i++) {
		g.addNode();
	}
	g.addEdge(0,1,0);
	g.addEdge(1,2,0);
	g.addEdge(2,4,0);
	g.addEdge(4,8,0);
	g.addEdge(8,16,0);
	g.addEdge(16,19,0);

	g.addEdge(3,5,0);
	g.addEdge(5,6,0);
	g.addEdge(6,7,0);
	g.addEdge(7,9,0);

	g.addEdge(10,11,0);
	g.addEdge(10,18,0);
	g.addEdge(10,12,0);
	g.addEdge(18,17,0);

	g.addEdge(13,14,0);

	// initialize ConnectedComponents
	ConnectedComponents ccs;
	ccs.run(g);

	// check result
	std::vector<node> comp = ccs.getComponent(0);
	EXPECT_TRUE(comp.size() == ccs.sizeOfComponent(0));
	EXPECT_TRUE(comp.at(0) == 0);
	EXPECT_TRUE(comp.at(1) == 1);
	EXPECT_TRUE(comp.at(2) == 2);
	EXPECT_TRUE(comp.at(3) == 4);
	EXPECT_TRUE(comp.at(4) == 8);
	EXPECT_TRUE(comp.at(5) == 16);
	EXPECT_TRUE(comp.at(6) == 19);
}

} /* namespace NetworKit */

#endif /*NOGTEST */

