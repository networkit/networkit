#include "GraphToolsGTest.h"
#include "../Graph.h"
#include "../GraphTools.h"

namespace NetworKit {

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

	auto Gcompact = GraphTools::getCompactedGraph(G);
	
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

	auto Gcompact = GraphTools::getCompactedGraph(G);
	
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

	auto Gcompact = GraphTools::getCompactedGraph(G);
	
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

	auto Gcompact = GraphTools::getCompactedGraph(G);
	
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

	auto Gcompact = GraphTools::getCompactedGraph(G);
	
	EXPECT_EQ(G.totalEdgeWeight(),Gcompact.totalEdgeWeight());
	EXPECT_NE(G.upperNodeIdBound(),Gcompact.upperNodeIdBound());
	EXPECT_EQ(G.numberOfNodes(),Gcompact.numberOfNodes());
	EXPECT_EQ(G.numberOfEdges(),Gcompact.numberOfEdges());
	EXPECT_EQ(G.isDirected(),Gcompact.isDirected());
	EXPECT_EQ(G.isWeighted(),Gcompact.isWeighted());
	// TODOish: find a deeper test to check if the structure of the graphs are the same, 
	// probably compare results of some algorithms or compare each edge with a reference node id map.
}

}