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

	Ghouse = bHouse.toGraph(useParallel());
}

TEST_P(GraphBuilderGTest, testEmptyGraph) {
	GraphBuilder b = createGraphBuilder();
	ASSERT_EQ(0u, b.numberOfNodes());

	Graph G = b.toGraph(useParallel());

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

	Graph G = b.toGraph(useParallel());

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
	if (isDirected()) {
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

TEST_P(GraphBuilderGTest, testDegreeIn) {
	if (isDirected()) {
		ASSERT_EQ(1u, this->Ghouse.degreeIn(0));
		ASSERT_EQ(2u, this->Ghouse.degreeIn(1));
		ASSERT_EQ(2u, this->Ghouse.degreeIn(2));
		ASSERT_EQ(1u, this->Ghouse.degreeIn(3));
		ASSERT_EQ(2u, this->Ghouse.degreeIn(4));
	} else {
		ASSERT_EQ(2u, this->Ghouse.degreeIn(0));
		ASSERT_EQ(4u, this->Ghouse.degreeIn(1));
		ASSERT_EQ(4u, this->Ghouse.degreeIn(2));
		ASSERT_EQ(3u, this->Ghouse.degreeIn(3));
		ASSERT_EQ(3u, this->Ghouse.degreeIn(4));
	}
}

TEST_P(GraphBuilderGTest, testDegreeOut) {
	if (isDirected()) {
		ASSERT_EQ(1u, this->Ghouse.degreeOut(0));
		ASSERT_EQ(2u, this->Ghouse.degreeOut(1));
		ASSERT_EQ(2u, this->Ghouse.degreeOut(2));
		ASSERT_EQ(2u, this->Ghouse.degreeOut(3));
		ASSERT_EQ(1u, this->Ghouse.degreeOut(4));
	} else {
		ASSERT_EQ(2u, this->Ghouse.degreeOut(0));
		ASSERT_EQ(4u, this->Ghouse.degreeOut(1));
		ASSERT_EQ(4u, this->Ghouse.degreeOut(2));
		ASSERT_EQ(3u, this->Ghouse.degreeOut(3));
		ASSERT_EQ(3u, this->Ghouse.degreeOut(4));
	}
}


/** EDGE MODIFIERS **/

TEST_P(GraphBuilderGTest, testAddEdge) {
	GraphBuilder b = createGraphBuilder(3);

	// Graph with 2 normal edges
	b.addEdge(0, 1, 4.51);
	b.addEdge(1, 2, 2.39);

	Graph G = b.toGraph(useParallel());

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


// /** GLOBAL PROPERTIES **/

TEST_P(GraphBuilderGTest, testIsWeighted) {
	ASSERT_EQ(isWeighted(), this->Ghouse.isWeighted());
}

TEST_P(GraphBuilderGTest, testIsDirected) {
	ASSERT_EQ(isDirected(), this->Ghouse.isDirected());
}

TEST_P(GraphBuilderGTest, testNumberOfSelfLoops) {
	GraphBuilder b = createGraphBuilder(3);
	b.addEdge(0, 1);
	b.addEdge(0, 0);
	Graph G = b.toGraph(useParallel());
	ASSERT_EQ(1u, G.numberOfSelfLoops());
}

TEST_P(GraphBuilderGTest, testUpperNodeIdBound) {
	ASSERT_EQ(5u, this->Ghouse.upperNodeIdBound());
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

		Graph G = b.toGraph(useParallel());
		
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


// /** Collections **/

// TEST_P(GraphBuilderGTest, testNodes) {
// 	GraphBuilder b = createGraphBuilder(3);
// 	G.addNode();
// 	G.removeNode(2);
// 	G.addNode();
// 	auto nodes = G.nodes();

// 	auto containsNode = [&nodes](node v) {
// 		return std::find(nodes.begin(), nodes.end(), v) != nodes.end();
// 	};

// 	ASSERT_EQ(G.numberOfNodes(), nodes.size());
// 	for (node v : nodes) {
// 		ASSERT_TRUE(containsNode(v));
// 	}
// }

// TEST_P(GraphBuilderGTest, testNeighbors) {
// 	auto neighbors = this->Ghouse.neighbors(1);
// 	auto containsNode = [&neighbors](node v) {
// 		return std::find(neighbors.begin(), neighbors.end(), v) != neighbors.end();
// 	};
// 	if(this->Ghouse.isDirected()){	
// 		ASSERT_TRUE(containsNode(0));
// 		ASSERT_TRUE(containsNode(4));
// 	}else{
// 		ASSERT_TRUE(containsNode(0));
// 		ASSERT_TRUE(containsNode(2));
// 		ASSERT_TRUE(containsNode(3));
// 		ASSERT_TRUE(containsNode(4));
// 	}

// }

// TEST_P(GraphBuilderGTest, testEdges) {
// 	// add self-loop
// 	this->Ghouse.addEdge(3, 3);
// 	auto isCorrectEdge = [&](node u, node v) {
// 		if (u == 3 && v == 3) {
// 			return true;
// 		}
// 		auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), std::make_pair(u, v));
// 		if (it != this->houseEdgesOut.end()) {
// 			return true;
// 		} else if (!this->Ghouse.isDirected()) {
// 			it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), std::make_pair(v, u));
// 			return it != this->houseEdgesOut.end();
// 		}
// 		return false;
// 	};
	
// 	auto edges = this->Ghouse.edges();
// 	ASSERT_EQ(this->m_house + 1, edges.size()); // plus self-loop
// 	for (auto e : edges) {
// 		ASSERT_TRUE(isCorrectEdge(e.first, e.second)) << "(" << e.first << ", " << e.second << ") is in edge array, but is not an edge of Ghouse";
// 	}
// }

// /** NODE ITERATORS **/

// TEST_P(GraphBuilderGTest, testForNodes) {
// 	GraphBuilder b = createGraphBuilder(3);
// 	std::vector<bool> visited(4, false);
// 	G.forNodes([&](node v) {
// 		ASSERT_FALSE(visited[v]);
// 		if (v == 2) {
// 			G.addNode();
// 		}
// 		visited[v] = true;
// 	});
// 	for (bool b : visited) {
// 		ASSERT_TRUE(b);
// 	}
// }

// TEST_P(GraphBuilderGTest, testParallelForNodes) {
// 	std::vector<node> visited;
// 	this->Ghouse.parallelForNodes([&](node u){
// 		#pragma omp critical
// 		visited.push_back(u);
// 	});
	
// 	std::sort(visited.begin(), visited.end());

// 	ASSERT_EQ(5u, visited.size());
// 	for (index i = 0; i < this->Ghouse.upperNodeIdBound(); i++) {
// 		ASSERT_EQ(i, visited[i]);
// 	}
// }






} /* namespace NetworKit */

#endif /*NOGTEST */
