/*
 * BasicGraph.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include "BasicGraphGTest.h"
#include "../BasicGraph.h"

namespace NetworKit {

using testing::Types;
typedef Types< Graph, WeightedGraph, DirectedGraph, WeightedDirectedGraph > graphImplementations;

TYPED_TEST_CASE(BasicGraphGTest, graphImplementations);

template <typename T>
void BasicGraphGTest<T>::SetUp() {
	/*
	 *    0
	 *   . \
	 *  /   \
	 * /     .
	 * 1 <-- 2
	 * ^ \  .|
	 * |  \/ |
	 * | / \ |
	 * |/   ..
	 * 3 <-- 4
	 *
	 * move you pen from node to node:
	 * 3 -> 1 -> 0 -> 2 -> 1 -> 4 -> 3 -> 2 -> 4
	 */
	Ghouse = T(5);
	houseEdgesOut = {
		{0, 2},
		{1, 0},
		{1, 4},
		{2, 1},
		{2, 4},
		{3, 1},
		{3, 2},
		{4, 3}
	};
	for (auto& e : houseEdgesOut) {
		Ghouse.addEdge(e.first, e.second);
	}
	n_house = 5;
	m_house = 8;
}

TYPED_TEST(BasicGraphGTest, getId) {
	TypeParam G1;
	TypeParam G2(5);

	ASSERT_TRUE(G1.getId() > 0);
	ASSERT_TRUE(G2.getId() > 0);	
	ASSERT_TRUE(G1.getId() < G2.getId());
}

TYPED_TEST(BasicGraphGTest, setName) {
	TypeParam G1(0);
	TypeParam G2(0);
	
	std::string s1 = "Graph 1";
	std::string s2 = "Graph 2";
	G1.setName(s1);
	G2.setName(s2);
	ASSERT_EQ(s1, G1.getName());
	ASSERT_EQ(s2, G2.getName());
}

TYPED_TEST(BasicGraphGTest, toString) {
	TypeParam G1(0);
	TypeParam G2(0);

	ASSERT_TRUE(G1.toString() != "");
	ASSERT_TRUE(G2.toString() != "");
}

TYPED_TEST(BasicGraphGTest, addNode) {
	TypeParam G(0);

	ASSERT_FALSE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(0u, G.numberOfNodes());

	G.addNode();
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(1u, G.numberOfNodes());

	TypeParam G2(2);
	ASSERT_TRUE(G2.hasNode(0));
	ASSERT_TRUE(G2.hasNode(1));
	ASSERT_FALSE(G2.hasNode(2));
	ASSERT_EQ(2u, G2.numberOfNodes());

	G2.addNode();
	G2.addNode();
	ASSERT_TRUE(G2.hasNode(2));
	ASSERT_TRUE(G2.hasNode(3));
	ASSERT_FALSE(G2.hasNode(4));
	ASSERT_EQ(4u, G2.numberOfNodes());
}

TYPED_TEST(BasicGraphGTest, removeNode) {
	TypeParam G(4);
	G.addEdge(0, 1);

	ASSERT_EQ(4u, G.numberOfNodes());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	EXPECT_ANY_THROW(G.removeNode(0));
	EXPECT_ANY_THROW(G.removeNode(1));

	G.removeNode(2);

	ASSERT_EQ(3u, G.numberOfNodes());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	G.removeNode(3);

	ASSERT_EQ(2u, G.numberOfNodes());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_FALSE(G.hasNode(3));
}

TYPED_TEST(BasicGraphGTest, addEdge) {
	TypeParam G(3);

	ASSERT_EQ(0u, G.numberOfEdges());
	ASSERT_FALSE(G.hasEdge(0, 0));
	ASSERT_FALSE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	G.addEdge(0, 1);

	ASSERT_EQ(1u, G.numberOfEdges());
	ASSERT_FALSE(G.hasEdge(0, 0));
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	G.addEdge(0, 0);

	ASSERT_EQ(2u, G.numberOfEdges());
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	// TODO weights?
}

TYPED_TEST(BasicGraphGTest, removeEdge) {
	TypeParam G(3);

	G.addEdge(0, 1);
	G.addEdge(0, 0);

	ASSERT_EQ(2u, G.numberOfEdges());
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	G.removeEdge(0, 1);
	ASSERT_EQ(1u, G.numberOfEdges());
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_FALSE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	// TODO weights?
}

TYPED_TEST(BasicGraphGTest, numberOfNodes) {
	ASSERT_EQ(this->n_house, this->Ghouse.numberOfNodes());

	TypeParam G1(0);
	ASSERT_EQ(0u, G1.numberOfNodes());
	G1.addNode();
	ASSERT_EQ(1u, G1.numberOfNodes());
	G1.addNode();
	ASSERT_EQ(2u, G1.numberOfNodes());
	G1.removeNode(0);
	ASSERT_EQ(1u, G1.numberOfNodes());
	G1.removeNode(1);
	ASSERT_EQ(0u, G1.numberOfNodes());
}

TYPED_TEST(BasicGraphGTest, numberOfEdges) {
	ASSERT_EQ(this->m_house, this->Ghouse.numberOfEdges());

	TypeParam G1(5);
	ASSERT_EQ(0u, G1.numberOfEdges());
	G1.addEdge(0, 1);
	ASSERT_EQ(1u, G1.numberOfEdges());
	G1.addEdge(1, 2);
	ASSERT_EQ(2u, G1.numberOfEdges());
	G1.removeEdge(0, 1);
	ASSERT_EQ(1u, G1.numberOfEdges());
	G1.removeEdge(1, 2);
	ASSERT_EQ(0u, G1.numberOfEdges());
}

TYPED_TEST(BasicGraphGTest, upperNodeIdBound) {
	ASSERT_EQ(5u, this->Ghouse.upperNodeIdBound());

	TypeParam G1(0);
	ASSERT_EQ(0u, G1.upperNodeIdBound());
	G1.addNode();
	ASSERT_EQ(1u, G1.upperNodeIdBound());
	G1.addNode();
	ASSERT_EQ(2u, G1.upperNodeIdBound());
	G1.removeNode(1);
	ASSERT_EQ(2u, G1.upperNodeIdBound());
	G1.addNode();
	ASSERT_EQ(3u, G1.upperNodeIdBound());
}

TYPED_TEST(BasicGraphGTest, degree) {
	if (this->Ghouse.isDirected()) {
		ASSERT_EQ(1u, this->Ghouse.degree(0));
		ASSERT_EQ(2u, this->Ghouse.degree(1));
		ASSERT_EQ(2u, this->Ghouse.degree(2));
		ASSERT_EQ(2u, this->Ghouse.degree(3));
		ASSERT_EQ(1u, this->Ghouse.degree(4));
	} else {
		ASSERT_EQ(2u, this->Ghouse.degree(0));
		ASSERT_EQ(4u, this->Ghouse.degree(1));
		ASSERT_EQ(4u, this->Ghouse.degree(2));
		ASSERT_EQ(3u, this->Ghouse.degree(3));
		ASSERT_EQ(3u, this->Ghouse.degree(4));
	}
}

TYPED_TEST(BasicGraphGTest, isEmpty) {
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

TYPED_TEST(BasicGraphGTest, weight) {
	if (this->Ghouse.isWeighted()) {
		// TODO
	} else {
		count c = 0;
		for (node u = 0; u < this->Ghouse.upperNodeIdBound(); u++) {
			for (node v = 0; v < this->Ghouse.upperNodeIdBound(); v++) {
				if (this->Ghouse.hasEdge(u, v)) {
					c++;
					ASSERT_EQ(defaultEdgeWeight, this->Ghouse.weight(u, v));
				} else {
					ASSERT_EQ(nullWeight, this->Ghouse.weight(u, v));
				}
			}
		}

		if (this->Ghouse.isDirected()) {
			ASSERT_EQ(this->m_house, c);
		} else {
			ASSERT_EQ(2 * this->m_house, c);
		}

		auto& e1 = this->houseEdgesOut[3];
		auto& e2 = this->houseEdgesOut[0];
		this->Ghouse.setWeight(e1.first, e1.second, 3.1415);
		this->Ghouse.setWeight(e2.first, e2.second, 2.71828);

		for (node u = 0; u < this->Ghouse.upperNodeIdBound(); u++) {
			for (node v = 0; v < this->Ghouse.upperNodeIdBound(); v++) {
				if (this->Ghouse.hasEdge(u, v)) {
					c++;
					ASSERT_EQ(defaultEdgeWeight, this->Ghouse.weight(u, v));
				} else {
					ASSERT_EQ(nullWeight, this->Ghouse.weight(u, v));
				}
			}
		}

		if (this->Ghouse.isDirected()) {
			ASSERT_EQ(2 * this->m_house, c);
		} else {
			ASSERT_EQ(4 * this->m_house, c);
		}
	}
}

TYPED_TEST(BasicGraphGTest, weightedDegree) {
	if (this->Ghouse.isWeighted()) {
		// TODO
	} else {
		if (this->Ghouse.isDirected()) {
			ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
			ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
		} else {
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
			ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
			ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
			ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
			ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
		}
	}
}

TYPED_TEST(BasicGraphGTest, totalEdgeWeight) {
	if (this->Ghouse.isWeighted()) {
		// TODO
	} else {
		ASSERT_EQ(this->m_house * defaultEdgeWeight, this->Ghouse.totalEdgeWeight());
	}
}

TYPED_TEST(BasicGraphGTest, forNodes) {
	TypeParam G(3);
	std::vector<bool> visited(4, false);
	G.forNodes([&](node v) {
		ASSERT_FALSE(visited[v]);
		if (v == 2) {
			G.addNode();
		}
		visited[v] = true;
	});
	for (bool b : visited) {
		ASSERT_TRUE(b);
	}
}

} /* namespace NetworKit */

#endif /*NOGTEST */
