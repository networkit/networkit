/*
 * GraphToolsGTest.cpp
 *
 *  Created on: 22.11.14
 *      Author: Maximilian Vogel
 */

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class GraphToolsGTest: public testing::Test {};

TEST_F(GraphToolsGTest, testGetContinuousOnContinuous) {
	Graph G(10);
	auto nodeIds = GraphTools::getContinuousNodeIds(G);
	std::unordered_map<node,node> reference = {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{9,9}};
	EXPECT_EQ(reference,nodeIds);
}

TEST_F(GraphToolsGTest, testGetContinuousOnDeletedNodes1) {
	Graph G(10);
	G.removeNode(0);
	G.removeNode(1);
	G.removeNode(2);
	G.removeNode(3);
	G.removeNode(4);
	auto nodeIds = GraphTools::getContinuousNodeIds(G);
	std::unordered_map<node,node> reference = {{5,0},{6,1},{7,2},{8,3},{9,4}};
	EXPECT_EQ(reference,nodeIds);
}

TEST_F(GraphToolsGTest, testGetContinuousOnDeletedNodes2) {
	Graph G(10);
	G.removeNode(0);
	G.removeNode(2);
	G.removeNode(4);
	G.removeNode(6);
	G.removeNode(8);
	auto nodeIds = GraphTools::getContinuousNodeIds(G);
	std::unordered_map<node,node> reference = {{1,0},{3,1},{5,2},{7,3},{9,4}};
	EXPECT_EQ(reference,nodeIds);
}

TEST_F(GraphToolsGTest, testGetCompactedGraphUndirectedUnweighted1) {
	Graph G(10,false,false);
	G.addEdge(0,1);
	G.addEdge(2,1);
	G.addEdge(0,3);
	G.addEdge(2,4);
	G.addEdge(3,6);
	G.addEdge(4,8);
	G.addEdge(5,9);
	G.addEdge(3,7);
	G.addEdge(5,7);

	auto nodeMap = GraphTools::getContinuousNodeIds(G);
	auto Gcompact = GraphTools::getCompactedGraph(G,nodeMap);

	EXPECT_EQ(G.numberOfNodes(),Gcompact.numberOfNodes());
	EXPECT_EQ(G.numberOfEdges(),Gcompact.numberOfEdges());
	EXPECT_EQ(G.isDirected(),Gcompact.isDirected());
	EXPECT_EQ(G.isWeighted(),Gcompact.isWeighted());
	// TODOish: find a deeper test to check if the structure of the graphs are the same,
	// probably compare results of some algorithms or compare each edge with a reference node id map.
}

TEST_F(GraphToolsGTest, testGetCompactedGraphUndirectedUnweighted2) {
	Graph G(10,false,false);
	G.removeNode(0);
	G.removeNode(2);
	G.removeNode(4);
	G.removeNode(6);
	G.removeNode(8);
	G.addEdge(1,3);
	G.addEdge(5,3);
	G.addEdge(7,5);
	G.addEdge(7,9);
	G.addEdge(1,9);

	auto nodeMap = GraphTools::getContinuousNodeIds(G);
	auto Gcompact = GraphTools::getCompactedGraph(G,nodeMap);

	EXPECT_NE(G.upperNodeIdBound(),Gcompact.upperNodeIdBound());
	EXPECT_EQ(G.numberOfNodes(),Gcompact.numberOfNodes());
	EXPECT_EQ(G.numberOfEdges(),Gcompact.numberOfEdges());
	EXPECT_EQ(G.isDirected(),Gcompact.isDirected());
	EXPECT_EQ(G.isWeighted(),Gcompact.isWeighted());
	// TODOish: find a deeper test to check if the structure of the graphs are the same,
	// probably compare results of some algorithms or compare each edge with a reference node id map.
}

TEST_F(GraphToolsGTest, testGetCompactedGraphUndirectedWeighted1) {
	Graph G(10,true,false);
	G.removeNode(0);
	G.removeNode(2);
	G.removeNode(4);
	G.removeNode(6);
	G.removeNode(8);
	G.addEdge(1,3,0.2);
	G.addEdge(5,3,2132.351);
	G.addEdge(7,5,3.14);
	G.addEdge(7,9,2.7);
	G.addEdge(1,9,0.12345);

	auto nodeMap = GraphTools::getContinuousNodeIds(G);
	auto Gcompact = GraphTools::getCompactedGraph(G,nodeMap);

	EXPECT_EQ(G.totalEdgeWeight(),Gcompact.totalEdgeWeight());
	EXPECT_NE(G.upperNodeIdBound(),Gcompact.upperNodeIdBound());
	EXPECT_EQ(G.numberOfNodes(),Gcompact.numberOfNodes());
	EXPECT_EQ(G.numberOfEdges(),Gcompact.numberOfEdges());
	EXPECT_EQ(G.isDirected(),Gcompact.isDirected());
	EXPECT_EQ(G.isWeighted(),Gcompact.isWeighted());
	// TODOish: find a deeper test to check if the structure of the graphs are the same,
	// probably compare results of some algorithms or compare each edge with a reference node id map.
}

TEST_F(GraphToolsGTest, testGetCompactedGraphDirectedWeighted1) {
	Graph G(10,true,true);
	G.removeNode(0);
	G.removeNode(2);
	G.removeNode(4);
	G.removeNode(6);
	G.removeNode(8);
	G.addEdge(1,3,0.2);
	G.addEdge(5,3,2132.351);
	G.addEdge(7,5,3.14);
	G.addEdge(7,9,2.7);
	G.addEdge(1,9,0.12345);

	auto nodeMap = GraphTools::getContinuousNodeIds(G);
	auto Gcompact = GraphTools::getCompactedGraph(G,nodeMap);

	EXPECT_EQ(G.totalEdgeWeight(),Gcompact.totalEdgeWeight());
	EXPECT_NE(G.upperNodeIdBound(),Gcompact.upperNodeIdBound());
	EXPECT_EQ(G.numberOfNodes(),Gcompact.numberOfNodes());
	EXPECT_EQ(G.numberOfEdges(),Gcompact.numberOfEdges());
	EXPECT_EQ(G.isDirected(),Gcompact.isDirected());
	EXPECT_EQ(G.isWeighted(),Gcompact.isWeighted());
	// TODOish: find a deeper test to check if the structure of the graphs are the same,
	// probably compare results of some algorithms or compare each edge with a reference node id map.
}

TEST_F(GraphToolsGTest, testGetCompactedGraphDirectedUnweighted1) {
	Graph G(10,false,true);
	G.removeNode(0);
	G.removeNode(2);
	G.removeNode(4);
	G.removeNode(6);
	G.removeNode(8);
	G.addEdge(1,3);
	G.addEdge(5,3);
	G.addEdge(7,5);
	G.addEdge(7,9);
	G.addEdge(1,9);
	auto nodeMap = GraphTools::getContinuousNodeIds(G);
	auto Gcompact = GraphTools::getCompactedGraph(G,nodeMap);

	EXPECT_EQ(G.totalEdgeWeight(),Gcompact.totalEdgeWeight());
	EXPECT_NE(G.upperNodeIdBound(),Gcompact.upperNodeIdBound());
	EXPECT_EQ(G.numberOfNodes(),Gcompact.numberOfNodes());
	EXPECT_EQ(G.numberOfEdges(),Gcompact.numberOfEdges());
	EXPECT_EQ(G.isDirected(),Gcompact.isDirected());
	EXPECT_EQ(G.isWeighted(),Gcompact.isWeighted());
	// TODOish: find a deeper test to check if the structure of the graphs are the same,
	// probably compare results of some algorithms or compare each edge with a reference node id map.
}

TEST_F(GraphToolsGTest, testInvertedMapping) {
	Graph G(10,false,true);
	G.removeNode(0);
	G.removeNode(2);
	G.removeNode(4);
	G.removeNode(6);
	G.removeNode(8);
	G.addEdge(1,3);
	G.addEdge(5,3);
	G.addEdge(7,5);
	G.addEdge(7,9);
	G.addEdge(1,9);
	auto nodeMap = GraphTools::getContinuousNodeIds(G);
	auto invertedNodeMap = GraphTools::invertContinuousNodeIds(nodeMap,G);

	EXPECT_EQ(6,invertedNodeMap.size());

	std::vector<node> reference = {1,3,5,7,9,10};
	EXPECT_EQ(reference,invertedNodeMap);
}

TEST_F(GraphToolsGTest, testRestoreGraph) {
	Graph G(10,false,true);
	G.removeNode(0);
	G.removeNode(2);
	G.removeNode(4);
	G.removeNode(6);
	G.removeNode(8);
	G.addEdge(1,3);
	G.addEdge(5,3);
	G.addEdge(7,5);
	G.addEdge(7,9);
	G.addEdge(1,9);
	auto nodeMap = GraphTools::getContinuousNodeIds(G);
	auto invertedNodeMap = GraphTools::invertContinuousNodeIds(nodeMap,G);
	std::vector<node> reference = {1,3,5,7,9,10};


	EXPECT_EQ(6,invertedNodeMap.size());
	EXPECT_EQ(reference,invertedNodeMap);

	auto Gcompact = GraphTools::getCompactedGraph(G,nodeMap);
	Graph Goriginal = GraphTools::restoreGraph(invertedNodeMap,Gcompact);

	EXPECT_EQ(Goriginal.totalEdgeWeight(),Gcompact.totalEdgeWeight());
	EXPECT_NE(Goriginal.upperNodeIdBound(),Gcompact.upperNodeIdBound());
	EXPECT_EQ(Goriginal.numberOfNodes(),Gcompact.numberOfNodes());
	EXPECT_EQ(Goriginal.numberOfEdges(),Gcompact.numberOfEdges());
	EXPECT_EQ(Goriginal.isDirected(),Gcompact.isDirected());
	EXPECT_EQ(Goriginal.isWeighted(),Gcompact.isWeighted());
}

TEST_F(GraphToolsGTest, testGetRemappedGraph) {
	for(bool directed : {false, true}) {
		const auto n = 4;
		Graph G(n, true, directed);
		for (auto i : {0, 1, 2})
			G.addEdge(i, i + 1, i);

		if (directed)
			G.addEdge(1, 1, 12);

		std::vector<node> perm(n);
		for (int i = 0; i < n; ++i) perm[i] = i;

		std::mt19937_64 gen;
		for (int iter = 0; iter < 10; iter++) {
			std::shuffle(perm.begin(), perm.end(), gen);
			auto G1 = GraphTools::getRemappedGraph(G, n, [&](node i) { return perm[i]; });
			ASSERT_EQ(G1.numberOfNodes(), n);
			ASSERT_EQ(G1.numberOfEdges(), G.numberOfEdges());
			ASSERT_EQ(G1.numberOfSelfLoops(), G.numberOfSelfLoops());

			for (int i = 0; i < n; ++i) {
				for (int j = 0; i < n; ++i) {
					ASSERT_EQ(G.hasEdge(i, j), G1.hasEdge(perm[i], perm[j]));
					ASSERT_EQ(G.weight(i, j), G1.weight(perm[i], perm[j]));
				}
			}
		}
	}
}

TEST_F(GraphToolsGTest, testGetRemappedGraphWithDelete) {
	for(bool directed : {false, true}) {
		const auto n = 4;
		Graph G(n, true, directed);
		for (auto i : {0, 1, 2})
			G.addEdge(i, i + 1, i);

		if (directed)
			G.addEdge(1, 1, 12);

		std::vector<node> perm(n);
		for (int i = 0; i < n; ++i) perm[i] = i;

		std::mt19937_64 gen;
		std::uniform_int_distribution<node> distr(0, n-1);
		for (int iter = 0; iter < 10; iter++) {
			std::shuffle(perm.begin(), perm.end(), gen);

			const auto del = distr(gen);

			auto G1 = GraphTools::getRemappedGraph(G, n,
				[&](node i) { return perm[i]; },
				[&](node i) { return i == del; }
			);

			auto expected_num_edges = G.numberOfEdges();
			expected_num_edges -= G.degree(del);
			if (directed)
				expected_num_edges -= G.degreeIn(del);
			//do double count self-loops
			expected_num_edges += G.hasEdge(del, del);

			ASSERT_EQ(G1.numberOfNodes(), n);
			ASSERT_EQ(G1.numberOfEdges(), expected_num_edges) << " del=" << del;
			ASSERT_EQ(G1.numberOfSelfLoops(), G.numberOfSelfLoops() - G.hasEdge(del, del)) << " del=" << del;

			for (int i = 0; i < n; ++i) {
				for (int j = 0; i < n; ++i) {
					if (i == static_cast<int>(del) || j == static_cast<int>(del)) {
						ASSERT_FALSE(G1.hasEdge(perm[i], perm[j])) << "i=" << i << " j=" << j << " del=" << del;
					} else {
						ASSERT_EQ(G.hasEdge(i, j), G1.hasEdge(perm[i], perm[j]));
						ASSERT_EQ(G.weight(i, j), G1.weight(perm[i], perm[j]));
					}
				}
			}
		}
	}
}

} // namespace NetworKit
