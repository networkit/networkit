/*
* EdmondsKarpGTest.cpp
 *
 *  Created on: 13.06.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "EdmondsKarpGTest.h"

namespace NetworKit {

TEST_F(EdmondsKarpGTest, tryEdmondsKarpP1) {
	Graph G(7, false);
	G.addEdge(0,1);
	G.addEdge(0,2);
	G.addEdge(0,3);
	G.addEdge(1,2);
	G.addEdge(1,4);
	G.addEdge(2,3);
	G.addEdge(2,4);
	G.addEdge(3,4);
	G.addEdge(3,5);
	G.addEdge(4,6);
	G.addEdge(5,6);

	EdmondsKarp edKa;
	int id;
	std::vector<node> sourceSet;
	edgeweight maxFlow = edKa.run(G, 0, 6, sourceSet, id);
	EXPECT_EQ(2, maxFlow) << "max flow is not correct";
	EXPECT_EQ(1, G.attribute_double(4,6,id));
	EXPECT_EQ(1, G.attribute_double(5,6,id));

	EXPECT_TRUE(std::find(sourceSet.begin(), sourceSet.end(), 0) != sourceSet.end());
	EXPECT_TRUE(std::find(sourceSet.begin(), sourceSet.end(), 1) != sourceSet.end());
	EXPECT_TRUE(std::find(sourceSet.begin(), sourceSet.end(), 2) != sourceSet.end());
	EXPECT_TRUE(std::find(sourceSet.begin(), sourceSet.end(), 3) != sourceSet.end());
	EXPECT_TRUE(std::find(sourceSet.begin(), sourceSet.end(), 4) != sourceSet.end());

	EXPECT_TRUE(std::find(sourceSet.begin(), sourceSet.end(), 5) == sourceSet.end());
	EXPECT_TRUE(std::find(sourceSet.begin(), sourceSet.end(), 6) == sourceSet.end());
}

TEST_F(EdmondsKarpGTest, tryEdmondsKarpP2) {
	Graph G(6, true);
	G.addEdge(0,1, 5);
	G.addEdge(0,2, 15);
	G.addEdge(1,3, 5);
	G.addEdge(1,4, 5);
	G.addEdge(2,3, 5);
	G.addEdge(2, 4, 5);
	G.addEdge(3,5, 15);
	G.addEdge(4,5, 5);
	EdmondsKarp edKa;
	edgeweight maxFlow = edKa.run(G, 0, 5);
	EXPECT_EQ(15, maxFlow) << "max flow is not correct";
}

TEST_F(EdmondsKarpGTest, tryEdmondsKarpUnconnected) {
	Graph G(6, true);
	G.addEdge(0,1, 5);
	G.addEdge(0,2, 15);
	G.addEdge(1,2, 5);
	G.addEdge(3, 4, 5);
	G.addEdge(3,5, 15);
	G.addEdge(4,5, 5);

	EdmondsKarp edKa;
	edgeweight maxFlow = edKa.run(G, 0, 5);
	EXPECT_EQ(0, maxFlow) << "max flow is not correct";
}

} /* namespace NetworKit */
