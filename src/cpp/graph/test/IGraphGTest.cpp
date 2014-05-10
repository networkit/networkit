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

TYPED_TEST_CASE_P(IGraphGTest);

TYPED_TEST_P(IGraphGTest, addNode) {
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

TYPED_TEST_P(IGraphGTest, removeNode) {
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

REGISTER_TYPED_TEST_CASE_P(IGraphGTest,
	addNode,
	removeNode
);

INSTANTIATE_TYPED_TEST_CASE_P(IGraph_Graph, IGraphGTest, Graph);

INSTANTIATE_TYPED_TEST_CASE_P(IGraph_DirectedGraph, IGraphGTest, DirectedGraph);

} /* namespace NetworKit */

#endif /*NOGTEST */
