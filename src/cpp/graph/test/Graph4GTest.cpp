/*
 * BasicGraph.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include "Graph4GTest.h"

namespace NetworKit {

INSTANTIATE_TEST_CASE_P(InstantiationName, Graph4GTest, testing::Values(
						std::make_tuple(false, false),
						std::make_tuple(true, false),
						std::make_tuple(false, true),
						std::make_tuple(true, true)));

bool Graph4GTest::isWeightedParameterized() const {
	return std::get<0>(GetParam());
}
bool Graph4GTest::isDirectedParameterized() const {
	return std::get<1>(GetParam());
}

Graph Graph4GTest::createParameterizedGraph(count n) const {
	bool weighted, directed;
	std::tie(weighted, directed) = GetParam();
	Graph G(n, weighted, directed);
	return G;
}

void Graph4GTest::SetUp() {
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

	Ghouse = createParameterizedGraph(5);
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
	Ahouse = {n_house, std::vector<edgeweight>(n_house, 0.0)};
	edgeweight ew = 1.0;
	for (auto& e : houseEdgesOut) {
		node u = e.first;
		node v = e.second;
		Ghouse.addEdge(u, v, ew);
		
		Ahouse[u][v] = ew;
		if (!Ghouse.isDirected()) {
			Ahouse[v][u] = ew;
		}
		
		if (Ghouse.isWeighted()) {
			ew += 1.0;
		}
	}
}

TEST_P(Graph4GTest, getId) {
	Graph G1 = createParameterizedGraph();
	Graph G2 = createParameterizedGraph(5);

	ASSERT_TRUE(G1.getId() > 0);
	ASSERT_TRUE(G2.getId() > 0);	
	ASSERT_TRUE(G1.getId() < G2.getId());
}

TEST_P(Graph4GTest, setName) {
	Graph G1 = createParameterizedGraph(0);
	Graph G2 = createParameterizedGraph(0);
	
	std::string s1 = "Graph 1";
	std::string s2 = "Graph 2";
	G1.setName(s1);
	G2.setName(s2);
	ASSERT_EQ(s1, G1.getName());
	ASSERT_EQ(s2, G2.getName());
}

TEST_P(Graph4GTest, toString) {
	Graph G1 = createParameterizedGraph(0);
	Graph G2 = createParameterizedGraph(0);

	ASSERT_TRUE(G1.toString() != "");
	ASSERT_TRUE(G2.toString() != "");
}

TEST_P(Graph4GTest, addNode) {
	Graph G = createParameterizedGraph();

	ASSERT_FALSE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(0u, G.numberOfNodes());

	G.addNode();
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(1u, G.numberOfNodes());

	Graph G2 = createParameterizedGraph(2);
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

TEST_P(Graph4GTest, removeNode) {
	Graph G = createParameterizedGraph();
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

TEST_P(Graph4GTest, addEdge) {
	Graph G = createParameterizedGraph();

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

TEST_P(Graph4GTest, removeEdge) {
	Graph G = createParameterizedGraph(3);

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

TEST_P(Graph4GTest, numberOfNodes) {
	ASSERT_EQ(this->n_house, this->Ghouse.numberOfNodes());

	Graph G1 = createParameterizedGraph(0);
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

TEST_P(Graph4GTest, numberOfEdges) {
	ASSERT_EQ(this->m_house, this->Ghouse.numberOfEdges());

	Graph G1 = createParameterizedGraph(5);
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

TEST_P(Graph4GTest, upperNodeIdBound) {
	ASSERT_EQ(5u, this->Ghouse.upperNodeIdBound());

	Graph G1 = createParameterizedGraph(0);
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

TEST_P(Graph4GTest, degree) {
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

TEST_P(Graph4GTest, weightedDegree) {
	if (this->Ghouse.isWeighted()) {
		if (this->Ghouse.isDirected()) {
			// only sum weight of outgoing edges
			ASSERT_EQ(1.0, this->Ghouse.weightedDegree(0));
			ASSERT_EQ(5.0, this->Ghouse.weightedDegree(1));
			ASSERT_EQ(9.0, this->Ghouse.weightedDegree(2));
			ASSERT_EQ(13.0, this->Ghouse.weightedDegree(3));
			ASSERT_EQ(8.0, this->Ghouse.weightedDegree(4));
		} else {
			ASSERT_EQ(3.0, this->Ghouse.weightedDegree(0));
			ASSERT_EQ(15.0, this->Ghouse.weightedDegree(1));
			ASSERT_EQ(17.0, this->Ghouse.weightedDegree(2));
			ASSERT_EQ(21.0, this->Ghouse.weightedDegree(3));
			ASSERT_EQ(16.0, this->Ghouse.weightedDegree(4));
		}
	} else {
		if (this->Ghouse.isDirected()) {
			// only count outgoing edges
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

TEST_P(Graph4GTest, isWeighted) {
	ASSERT_EQ(isWeightedParameterized(), this->Ghouse.isWeighted());
}

TEST_P(Graph4GTest, isDirected) {
	ASSERT_EQ(isDirectedParameterized(), this->Ghouse.isDirected());
}

TEST_P(Graph4GTest, isEmpty) {
	Graph G1 = createParameterizedGraph(0);
	Graph G2 = createParameterizedGraph(2);

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

TEST_P(Graph4GTest, weight) {
	this->Ghouse.forNodes([&](node u) {
		this->Ghouse.forNodes([&](node v) {
			ASSERT_EQ(this->Ahouse[u][v], this->Ghouse.weight(u, v));
		});
	});
}

TEST_P(Graph4GTest, totalEdgeWeight) {
	Graph G1 = createParameterizedGraph(5);
	Graph G2 = createParameterizedGraph(5);
	G2.addEdge(0, 1, 3.14);

	if (this->Ghouse.isWeighted()) {
		ASSERT_EQ(0.0, G1.totalEdgeWeight());
		ASSERT_EQ(3.14, G2.totalEdgeWeight());
		ASSERT_EQ(36.0, this->Ghouse.totalEdgeWeight());
	} else {
		ASSERT_EQ(0 * defaultEdgeWeight, G1.totalEdgeWeight());
		ASSERT_EQ(1 * defaultEdgeWeight, G2.totalEdgeWeight());
		ASSERT_EQ(8 * defaultEdgeWeight, this->Ghouse.totalEdgeWeight());
	}
}

TEST_P(Graph4GTest, forNodes) {
	Graph G = createParameterizedGraph(3);
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

TEST_P(Graph4GTest, degreeInOut) {
	if (isDirectedParameterized()) {
		std::vector<count> inDegrees(this->Ghouse.upperNodeIdBound(), 0);
		std::vector<count> outDegrees(this->Ghouse.upperNodeIdBound(), 0);
		for (auto& e : this->houseEdgesOut) {
			outDegrees[e.first]++;
			inDegrees[e.second]++;
		}

		this->Ghouse.forNodes([&](node v) {
			EXPECT_EQ(inDegrees[v], this->Ghouse.degreeIn(v));
			EXPECT_EQ(outDegrees[v], this->Ghouse.degreeOut(v));
			EXPECT_EQ(outDegrees[v], this->Ghouse.degree(v));
		});
	} else {
		// TODO
	}
}

TEST_P(Graph4GTest, forEdgesOf) {
	if (isDirectedParameterized()) {
		count m = 0;
		std::vector<bool> visited(this->m_house, false);

		this->Ghouse.forNodes([&](node u) {
			this->Ghouse.forEdgesOf(u, [&](node v, node w) {
				// edges should be v to w, so if we iterate over edges from u, u should be equal v
				EXPECT_EQ(u, v);
				
				auto e = std::make_pair(v, w);
				// find edge
				auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);

				EXPECT_FALSE(it == this->houseEdgesOut.end()); // check if edge is allowed to exists
			
				// find index in edge array
				int i = std::distance(this->houseEdgesOut.begin(), it);
				EXPECT_FALSE(visited[i]); // make sure edge was not visited before (would be visited twice)
				
				// mark edge as visited
				visited[i] = true;
				m++;
			});
		});

		EXPECT_EQ(this->m_house, m);
		for (auto b : visited) {
			EXPECT_TRUE(b);
		}
	} else {
		// TODO
	}
}

TEST_P(Graph4GTest, BFSfrom) {
	if (isDirectedParameterized()) {
		std::vector<count> visitedOrder(5, none);
		index i = 0;
		this->Ghouse.BFSfrom(3, [&](node v) {
			EXPECT_EQ(none, visitedOrder[v]); // visit every node once
			visitedOrder[v] = i++;
		});
		// have we visited all nodes
		for (count l : visitedOrder) {
			EXPECT_TRUE(l != none);
		}

		// root on level 0
		EXPECT_EQ(0u, visitedOrder[3]);

		// level 1
		EXPECT_TRUE( (visitedOrder[1] == 1) ^ (visitedOrder[1] == 2) );
		EXPECT_TRUE( (visitedOrder[2] == 1) ^ (visitedOrder[2] == 2) );

		// level 2
		EXPECT_TRUE( (visitedOrder[0] == 3) ^ (visitedOrder[0] == 4) );
		EXPECT_TRUE( (visitedOrder[4] == 3) ^ (visitedOrder[4] == 4) );
	} else {
		// TODO
	}
}

TEST_P(Graph4GTest, DFSfrom) {
	if (isDirectedParameterized()) {
		std::vector<count> visitedOrder(5, none);
		index i = 0;
		this->Ghouse.DFSfrom(3, [&](node v) {
			EXPECT_EQ(none, visitedOrder[v]); // visit every node once
			visitedOrder[v] = i++;
		});

		// have we visited all nodes
		for (count l : visitedOrder) {
			EXPECT_TRUE(l != none);
		}

		// root on level 0
		EXPECT_EQ(0u, visitedOrder[3]);

		// level 1
		EXPECT_TRUE( (visitedOrder[1] == 1) ^ (visitedOrder[2] == 1) );

		// level 2
		EXPECT_TRUE( (visitedOrder[0] == 2) ^ (visitedOrder[1] == 2) ^ (visitedOrder[4] == 2) );

		// level 3
		EXPECT_TRUE( (visitedOrder[2] == 3) ^ (visitedOrder[0] == 3) ^ (visitedOrder[4] == 3) ^ (visitedOrder[1] == 3) );

		// level 4
		EXPECT_TRUE( (visitedOrder[2] == 4) ^ (visitedOrder[4] == 4) ^ (visitedOrder[0] == 4) );
	} else {
		// TODO
	}
}


} /* namespace NetworKit */

#endif /*NOGTEST */
