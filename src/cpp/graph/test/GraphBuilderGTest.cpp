/*
 * GraphBuilderGTest.cpp
 *
 *  Created on: 14.08.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include <algorithm>

#include "GraphBuilderGTest.h"

namespace NetworKit {

INSTANTIATE_TEST_CASE_P(InstantiationName, GraphBuilderGTest, testing::Values(
						std::make_tuple(false, false, false),
						std::make_tuple(true, false, false),
						std::make_tuple(false, true, false),
						std::make_tuple(true, true, false),
						std::make_tuple(false, false, true),
						std::make_tuple(true, false, true),
						std::make_tuple(false, true, true),
						std::make_tuple(true, true, true)
						));

bool GraphBuilderGTest::isWeighted() const {
	return std::get<0>(GetParam());
}
bool GraphBuilderGTest::isDirected() const {
	return std::get<1>(GetParam());
}

bool GraphBuilderGTest::useParallel() const {
	return std::get<2>(GetParam());
}

GraphBuilder GraphBuilderGTest::createGraphBuilder(count n) const {
	bool weighted, directed, parallel;
	std::tie(weighted, directed, parallel) = GetParam();
	GraphBuilder b(n, weighted, directed);
	return b;
}

Graph GraphBuilderGTest::toGraph(GraphBuilder& b) const { return b.toGraph(useParallel()); };

void GraphBuilderGTest::SetUp() {
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
	n_house = 5;
	m_house = 8;

	bHouse = createGraphBuilder(5);
	houseEdgesOut = {
		{3, 1},
		{1, 0},
		{0, 2},
		{2, 1},
		{1, 4},
		{4, 3},
		{3, 2},
		{2, 4}
	};
	Ahouse = {n_house, std::vector<edgeweight>(n_house, 0.0)};
	edgeweight ew = 1.0;
	for (auto& e : houseEdgesOut) {
		node u = e.first;
		node v = e.second;
		bHouse.addEdge(u, v, ew);
		
		Ahouse[u][v] = ew;
	
		if (!bHouse.isDirected()) {
			Ahouse[v][u] = ew;
		}
		
		if (bHouse.isWeighted()) {
			ew += 1.0;
		}
	}
}

TEST_P(GraphBuilderGTest, testEmptyGraph) {
	GraphBuilder b = createGraphBuilder();
	ASSERT_EQ(0u, b.numberOfNodes());

	Graph G = toGraph(b);

	ASSERT_EQ(0u, G.numberOfNodes());
	ASSERT_EQ(0u, G.numberOfEdges());
	ASSERT_TRUE(G.isEmpty());
}

TEST_P(GraphBuilderGTest, testAddNode) {
	GraphBuilder b = createGraphBuilder();

	b.addNode();
	ASSERT_EQ(1u, b.numberOfNodes());

	b.addNode();
	b.addNode();
	ASSERT_EQ(3u, b.numberOfNodes());

	Graph G = toGraph(b);

	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_FALSE(G.hasNode(3));
	ASSERT_EQ(3u, G.numberOfNodes());
	ASSERT_EQ(0u, G.numberOfEdges());
	ASSERT_FALSE(G.isEmpty());
}


/** NODE PROPERTIES **/

TEST_P(GraphBuilderGTest, testDegree) {
	Graph Ghouse = toGraph(this->bHouse);
	if (isDirected()) {
		ASSERT_EQ(1u, Ghouse.degree(0));
		ASSERT_EQ(2u, Ghouse.degree(1));
		ASSERT_EQ(2u, Ghouse.degree(2));
		ASSERT_EQ(2u, Ghouse.degree(3));
		ASSERT_EQ(1u, Ghouse.degree(4));
	} else {
		ASSERT_EQ(2u, Ghouse.degree(0));
		ASSERT_EQ(4u, Ghouse.degree(1));
		ASSERT_EQ(4u, Ghouse.degree(2));
		ASSERT_EQ(3u, Ghouse.degree(3));
		ASSERT_EQ(3u, Ghouse.degree(4));
	}
}

TEST_P(GraphBuilderGTest, testDegreeIn) {
	Graph Ghouse = toGraph(this->bHouse);
	if (isDirected()) {
		ASSERT_EQ(1u, Ghouse.degreeIn(0));
		ASSERT_EQ(2u, Ghouse.degreeIn(1));
		ASSERT_EQ(2u, Ghouse.degreeIn(2));
		ASSERT_EQ(1u, Ghouse.degreeIn(3));
		ASSERT_EQ(2u, Ghouse.degreeIn(4));
	} else {
		ASSERT_EQ(2u, Ghouse.degreeIn(0));
		ASSERT_EQ(4u, Ghouse.degreeIn(1));
		ASSERT_EQ(4u, Ghouse.degreeIn(2));
		ASSERT_EQ(3u, Ghouse.degreeIn(3));
		ASSERT_EQ(3u, Ghouse.degreeIn(4));
	}
}

TEST_P(GraphBuilderGTest, testDegreeOut) {
	Graph Ghouse = toGraph(this->bHouse);
	if (isDirected()) {
		ASSERT_EQ(1u, Ghouse.degreeOut(0));
		ASSERT_EQ(2u, Ghouse.degreeOut(1));
		ASSERT_EQ(2u, Ghouse.degreeOut(2));
		ASSERT_EQ(2u, Ghouse.degreeOut(3));
		ASSERT_EQ(1u, Ghouse.degreeOut(4));
	} else {
		ASSERT_EQ(2u, Ghouse.degreeOut(0));
		ASSERT_EQ(4u, Ghouse.degreeOut(1));
		ASSERT_EQ(4u, Ghouse.degreeOut(2));
		ASSERT_EQ(3u, Ghouse.degreeOut(3));
		ASSERT_EQ(3u, Ghouse.degreeOut(4));
	}
}


/** EDGE MODIFIERS **/

TEST_P(GraphBuilderGTest, testAddEdge) {
	GraphBuilder b = createGraphBuilder(3);

	// Graph with 2 normal edges
	b.addEdge(0, 1, 4.51);
	b.addEdge(1, 2, 2.39);

	Graph G = toGraph(b);

	ASSERT_EQ(2u, G.numberOfEdges());
	ASSERT_FALSE(G.hasEdge(0, 2)); // was never added
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_TRUE(G.hasEdge(1, 2));
	ASSERT_FALSE(G.hasEdge(2, 2)); // will be added later

	// check weights
	if (G.isWeighted()) {
		ASSERT_EQ(4.51, G.weight(0, 1));
		ASSERT_EQ(2.39, G.weight(1, 2));
	} else {
		ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
		ASSERT_EQ(defaultEdgeWeight, G.weight(1, 2));
	}

	if (G.isDirected()) {
		ASSERT_FALSE(G.hasEdge(1, 0));
		ASSERT_FALSE(G.hasEdge(2, 1));

		// add edge in the other direction
		// note: bidirectional edges are not supported, so both edges have different weights
		G.addEdge(2, 1, 6.23);
		ASSERT_TRUE(G.hasEdge(2, 1));
		if (G.isWeighted()) {
			ASSERT_EQ(2.39, G.weight(1, 2));
			ASSERT_EQ(6.23, G.weight(2, 1));
		} else {
			ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
		}
	} else {
		ASSERT_TRUE(G.hasEdge(1, 0));
		ASSERT_TRUE(G.hasEdge(2, 1));
		if (G.isWeighted()) {
			ASSERT_EQ(4.51, G.weight(1, 0));
			ASSERT_EQ(2.39, G.weight(2, 1));
		} else {
			ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
			ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
		}	
	}
}


/** GLOBAL PROPERTIES **/

TEST_P(GraphBuilderGTest, testIsWeighted) {
	ASSERT_EQ(isWeighted(), this->bHouse.isWeighted());
	Graph Ghouse = toGraph(this->bHouse);
	ASSERT_EQ(isWeighted(), Ghouse.isWeighted());
}

TEST_P(GraphBuilderGTest, testIsDirected) {
	ASSERT_EQ(isDirected(), this->bHouse.isDirected());
	Graph Ghouse = toGraph(this->bHouse);
	ASSERT_EQ(isDirected(), Ghouse.isDirected());
}

TEST_P(GraphBuilderGTest, testNumberOfSelfLoops) {
	GraphBuilder b = createGraphBuilder(3);
	b.addEdge(0, 1);
	b.addEdge(0, 0);
	Graph G = toGraph(b);
	ASSERT_EQ(1u, G.numberOfSelfLoops());
}

TEST_P(GraphBuilderGTest, testUpperNodeIdBound) {
	ASSERT_EQ(5u, this->bHouse.upperNodeIdBound());
	Graph Ghouse = toGraph(this->bHouse);
	ASSERT_EQ(5u, Ghouse.upperNodeIdBound());
}


/** EDGE ATTRIBUTES **/

TEST_P(GraphBuilderGTest, testSetWeight) {
	GraphBuilder b = createGraphBuilder(10);
	b.addEdge(0, 1);
	b.addEdge(1, 2);

	if (isWeighted()) {
		// edges should get weight defaultWeight on creation and setWeight should overwrite this
		b.setWeight(1, 2, 2.718);
		
		// setting an edge weight should create the edge if it doesn't exists
		b.setWeight(5, 6, 56.0);
		// directed graphs are not symmetric, undirected are
		b.setWeight(3, 4, 2.718);
		b.setWeight(4, 3, 5.243);
		
		// self-loop
		b.addEdge(8, 8, 2.5);
		b.setWeight(8, 8, 3.14);

		Graph G = toGraph(b);

		// edges should get weight defaultWeight on creation and setWeight should overwrite this
		ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
		ASSERT_EQ(2.718, G.weight(1, 2));
		if (isDirected()) {
			ASSERT_EQ(nullWeight, G.weight(1, 0));
			ASSERT_EQ(nullWeight, G.weight(2, 1));
		} else {
			// undirected graph is symmetric
			ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
			ASSERT_EQ(2.718, G.weight(2, 1));
		}

		// setting an edge weight should create the edge if it doesn't exists
		ASSERT_EQ(56.0, G.weight(5, 6));
		ASSERT_EQ(isDirected() ? nullWeight : 56.0, G.weight(6, 5));
		ASSERT_TRUE(G.hasEdge(5, 6));

		// directed graphs are not symmetric, undirected are
		if (isDirected()) {
			ASSERT_EQ(2.718, G.weight(3, 4));
			ASSERT_EQ(5.243, G.weight(4, 3));
		} else {
			// we have actually 2 edges in both direction. It is not defined which weight is returned.
			ASSERT_TRUE(G.weight(3, 4) == 2.718 || G.weight(3, 4) == 5.243);
			ASSERT_TRUE(G.weight(4, 3) == 2.718 || G.weight(4, 3) == 5.243);
		}
		
		// self-loop
		ASSERT_EQ(3.14, G.weight(8, 8));
	} else {
		EXPECT_ANY_THROW(b.setWeight(0, 1, 1.5));
	}
}


/** toGraph **/

TEST_P(GraphBuilderGTest, testSameAsGraph) {
	FAIL();
}

TEST_P(GraphBuilderGTest, testForValidStateAfterToGraph) {
	Graph Ghouse = toGraph(this->bHouse);

	ASSERT_TRUE(this->bHouse.isEmpty());
	ASSERT_EQ(0u, this->bHouse.numberOfNodes());
	ASSERT_EQ(0u, this->bHouse.upperNodeIdBound());
	ASSERT_EQ(isWeighted(), this->bHouse.isWeighted());
	ASSERT_EQ(isDirected(), this->bHouse.isDirected());
	this->bHouse.forNodes([&](node v) {
		FAIL();
	});

	Graph G1 = toGraph(this->bHouse);
	ASSERT_TRUE(G1.isEmpty());
	ASSERT_EQ(0u, G1.numberOfNodes());
	ASSERT_EQ(0u, G1.upperNodeIdBound());
	ASSERT_EQ(isWeighted(), G1.isWeighted());
	ASSERT_EQ(isDirected(), G1.isDirected());

	node v = this->bHouse.addNode();
	node u = this->bHouse.addNode();
	this->bHouse.addEdge(v, u, 0.25);

	Graph G2 = toGraph(this->bHouse);
	ASSERT_FALSE(G2.isEmpty());
	ASSERT_EQ(2u, G2.numberOfNodes());
	ASSERT_EQ(1u, G2.numberOfEdges());
	ASSERT_TRUE(G2.hasEdge(v, u));
	if (!isDirected()) {
		ASSERT_TRUE(G2.hasEdge(u, v));
	}
	if (isWeighted()) {
		ASSERT_EQ(0.25, G2.weight(v, u));
	} else {
		ASSERT_EQ(1.0, G2.weight(v, u));
	}
}

/** NODE ITERATORS **/

TEST_P(GraphBuilderGTest, testForNodes) {
	GraphBuilder b = createGraphBuilder(3);
	std::vector<bool> visited(4, false);
	b.forNodes([&](node v) {
		ASSERT_FALSE(visited[v]);
		if (v == 2) {
			b.addNode();
		}
		visited[v] = true;
	});
	for (bool x : visited) {
		ASSERT_TRUE(x);
	}
}

TEST_P(GraphBuilderGTest, testParallelForNodes) {
	std::vector<node> visited;
	this->bHouse.parallelForNodes([&](node u) {
		#pragma omp critical
		visited.push_back(u);
	});
	
	std::sort(visited.begin(), visited.end());

	ASSERT_EQ(5u, visited.size());
	for (index i = 0; i < this->bHouse.upperNodeIdBound(); i++) {
		ASSERT_EQ(i, visited[i]);
	}
}

TEST_P(GraphBuilderGTest, testForNodePairs) {
	FAIL();
}

TEST_P(GraphBuilderGTest, testParallelForNodePairs) {
	FAIL();
}

} /* namespace NetworKit */

#endif /*NOGTEST */
