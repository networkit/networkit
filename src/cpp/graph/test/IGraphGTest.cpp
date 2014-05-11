/*
 * IGraphGTest.cpp
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include "IGraphGTest.h"

#include "../Graph.h"
#include "../DirectedGraph.h"

namespace NetworKit {

template <typename T>
IGraphGTest<T>::IGraphGTest() {
	// TODO Auto-generated constructor stub

}

template <typename T>
IGraphGTest<T>::~IGraphGTest() {
	// TODO Auto-generated destructor stub
}

using testing::Types;
typedef Types<Graph, DirectedGraph> graphImplementations;

TYPED_TEST_CASE(IGraphGTest, graphImplementations);

TYPED_TEST(IGraphGTest, idNameAndToString) {
	TypeParam G1(0);
	TypeParam G2(0);
	
	ASSERT_TRUE(G1.getId() < G2.getId());
	
	std::string s1 = "Graph 1";
	std::string s2 = "Graph 2";
	G1.setName(s1);
	G2.setName(s2);
	ASSERT_EQ(G1.getName(), s1);
	ASSERT_EQ(G2.getName(), s2);

	ASSERT_TRUE(G1.toString() != "");
	ASSERT_TRUE(G2.toString() != "");
}

TYPED_TEST(IGraphGTest, addNode) {
	TypeParam G(0);

	ASSERT_FALSE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(G.numberOfNodes(), 0);

	G.addNode();
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(G.numberOfNodes(), 1);

	TypeParam G2(2);
	ASSERT_TRUE(G2.hasNode(0));
	ASSERT_TRUE(G2.hasNode(1));
	ASSERT_FALSE(G2.hasNode(2));
	ASSERT_EQ(G2.numberOfNodes(), 2);

	G2.addNode();
	G2.addNode();
	ASSERT_TRUE(G2.hasNode(2));
	ASSERT_TRUE(G2.hasNode(3));
	ASSERT_FALSE(G2.hasNode(4));
	ASSERT_EQ(G2.numberOfNodes(), 4);
}

TYPED_TEST(IGraphGTest, removeNode) {
	TypeParam G(4);
	G.addEdge(0, 1);

	ASSERT_EQ(G.numberOfNodes(), 4);
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	EXPECT_ANY_THROW(G.removeNode(0));
	EXPECT_ANY_THROW(G.removeNode(1));

	G.removeNode(2);

	ASSERT_EQ(G.numberOfNodes(), 3);
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	G.removeNode(3);

	ASSERT_EQ(G.numberOfNodes(), 2);
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_FALSE(G.hasNode(3));
}

TYPED_TEST(IGraphGTest, addEdge) {
	TypeParam G(3);

	ASSERT_EQ(G.numberOfEdges(), 0);
	ASSERT_FALSE(G.hasEdge(0, 0));
	ASSERT_FALSE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	G.addEdge(0, 1);

	ASSERT_EQ(G.numberOfEdges(), 1);
	ASSERT_FALSE(G.hasEdge(0, 0));
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	G.addEdge(0, 0);

	ASSERT_EQ(G.numberOfEdges(), 2);
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));
}

TYPED_TEST(IGraphGTest, removeEdge) {
	TypeParam G(3);

	G.addEdge(0, 1);
	G.addEdge(0, 0);

	ASSERT_EQ(G.numberOfEdges(), 2);
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	G.removeEdge(0, 1);
	ASSERT_EQ(G.numberOfEdges(), 1);
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_FALSE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));
}

TYPED_TEST(IGraphGTest, isEmpty) {
	TypeParam G1(0);
	TypeParam G2(2);

	ASSERT_TRUE(G1.isEmpty());
	ASSERT_FALSE(G2.isEmpty());

	node v = G1.addNode();
	G2.removeNode(G2.randomNode());
	ASSERT_FALSE(G1.isEmpty());
	ASSERT_FALSE(G2.isEmpty());

	G1.removeNode(v);
	G2.removeNode(G2.randomNode());
	ASSERT_TRUE(G1.isEmpty());
	ASSERT_TRUE(G2.isEmpty());
}

} /* namespace NetworKit */

#endif /*NOGTEST */
