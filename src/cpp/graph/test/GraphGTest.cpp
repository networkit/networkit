/*
 * GraphGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "GraphGTest.h"
#include "../BFS.h"
#include "../Dijkstra.h"

namespace NetworKit {

// TEST_F(GraphGTest, minMaxDegree) {
// 	Graph G = this->gen.makeRandomGraph(100, 0.1);

// 	node argmin = 0;
// 	node argmax = 0;
// 	count min = G.degree(argmin);
// 	count max = G.degree(argmax);

// 	G.forNodes([&](node v) {
// 		if (G.degree(v) < min) {
// 			argmin = v;
// 			min = G.degree(v);
// 		} else if (G.degree(v) > max) {
// 			argmax = v;
// 			max = G.degree(v);
// 		}
// 	});

// 	ASSERT_EQ(min, G.minDegree());
// 	ASSERT_EQ(argmin, G.argminDegree());
// 	ASSERT_EQ(max, G.maxDegree());
// 	ASSERT_EQ(argmax, G.argmaxDegree());
// }


//TEST_F(GraphGTest, testEdgeIteration) {
//
//	int64_t n = 100;
//	Graph G = this->gen.makeCompleteGraph(n);
//
//	int64_t edgeCount = 0;
//	READ_ONLY_FORALL_EDGES_BEGIN(G) {
//		node u = EDGE_SOURCE;
//		node v = EDGE_DEST;
//		if (u < v) {
//			edgeCount += 1;
//		}
//	} READ_ONLY_FORALL_EDGES_END();
//
//	EXPECT_EQ((n * (n-1)) / 2, edgeCount) << "There are (n * (n-1)) / 2 undirected edges in a compete graph";
//}


TEST_F(GraphGTest, testLambdaEdgeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t edgeCount = 0;
	G.forEdges([&](node u, node v) {
		if (u > v) {
			edgeCount += 1;
		}
	});

	EXPECT_EQ((n * (n-1)) / 2, edgeCount) << "There are (n * (n-1)) / 2 undirected edges in a compete graph";
}



TEST_F(GraphGTest, testParallelLambdaEdgeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t edgeCount = 0;
	G.parallelForEdges([&](node u, node v) {
		if (u > v) {
			#pragma omp atomic update
			edgeCount += 1;
		}
	});

	EXPECT_EQ((n * (n-1)) / 2, edgeCount) << "There are (n * (n-1)) / 2 undirected edges in a compete graph";
}




TEST_F(GraphGTest, testLambdaEdgeModification) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	G.forEdges([&](node u, node v){
		G.removeEdge(u, v);
	});

	EXPECT_EQ(0u, G.numberOfEdges()) << "all edges should have been deleted";
}



TEST_F(GraphGTest, testLambdaNodeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t nodeCount = 0;
	G.forNodes([&](node v) {
		nodeCount++;
	});

	EXPECT_EQ(n, nodeCount);
}


TEST_F(GraphGTest, testParallelLambdaNodeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t nodeCount = 0;
	G.parallelForNodes([&](node v) {
		#pragma omp atomic update
		nodeCount++;
	});

	EXPECT_EQ(n, nodeCount);
}


TEST_F(GraphGTest, testLambdaNeighborIteration) {

	Graph G(4);
	node v = 0;
	G.addEdge(v, 1);
	G.addEdge(v, 2);
	G.addEdge(v, 3);

	int neighborCount = 0;
	G.forNeighborsOf(v, [&](node w){
		neighborCount += 1;
	});

	EXPECT_EQ(3, neighborCount) << "node v has 3 neighbors";
}


TEST_F(GraphGTest, testLambdaIncidentEdgeIteration) {

	Graph G(4);
	node v = 0;
	G.addEdge(v, 1);
	G.addEdge(v, 2);
	G.addEdge(v, 3);

	int edgeCount = 0;
	G.forEdgesOf(v, [&](node v, node w){
		edgeCount += 1;
	});

	EXPECT_EQ(3, edgeCount) << "node v has 3 incident edges";
}


TEST_F(GraphGTest, testNodeBFS) {

	Graph G(4);
	node v = 0;
	G.addEdge(v, 1);
	G.addEdge(v, 2);
	G.addEdge(v, 3);

	std::vector<int> visited(4, 0);

	int nodeCount = 0;
	G.BFSfrom(v, [&](node w) {
		nodeCount += 1;
	});

	EXPECT_EQ(4, nodeCount) << "4 nodes should have been visited by BFS";
}

//TEST_F(GraphGTest, testEdgeBFS) {
//
//	Graph G(4);
//	node v = 0;
//	G.addEdge(v, 1);
//	G.addEdge(v, 2);
//	G.addEdge(v, 3);
//	G.addEdge(3, 2);
//
//	int edgeCount = 0;
//	G.breadthFirstEdgesFrom(v, [&](node x, node y) {
//		edgeCount += 1;
//	});
//
//	EXPECT_EQ(4, edgeCount) << "4 edges should have been visited by BFS";
//}

TEST_F(GraphGTest, testDegree) {
	Graph G(5);

	G.addEdge(1, 4);
	G.addEdge(3, 4);
	G.removeNode(2);

	EXPECT_EQ(0u, G.degree(0)) << "degree of node 0 should be 0";
	EXPECT_EQ(1u, G.degree(1)) << "degree of node 1 should be 1";
	// EXPECT_ANY_THROW(G.degree(2)) << "node 2 should not exist";
	EXPECT_EQ(1u, G.degree(3)) << "degree of node 3 should be 1";
	EXPECT_EQ(2u, G.degree(4)) << "degree of node 4 should be 2";
	// EXPECT_ANY_THROW(G.degree(5)) << "node 2 should not exist";
	// EXPECT_ANY_THROW(G.degree(6)) << "node 2 should not exist";
}


TEST_F(GraphGTest, testNodeIteration) {
	int64_t n = 42;
	Graph G(n);

	int64_t counter = 0;

	G.forNodes([&](node v) {
		counter += 1;
	});

	EXPECT_EQ(n, counter) << "all nodes should have been counted";
}


TEST_F(GraphGTest, testNumberOfEdges) {
	int64_t n = 3;
	Graph G(n);

	G.addEdge(0, 1);
	EXPECT_EQ(1u, G.numberOfEdges()) << "G should have 1 edge now";

	G.addEdge(0, 2);
	EXPECT_EQ(2u, G.numberOfEdges()) << "G should have 2 edges now";

}


TEST_F(GraphGTest, testHasEdge) {
	int64_t n = 3;
	Graph G(n);

	G.addEdge(0, 1);
	EXPECT_TRUE(G.hasEdge(0, 1)) << "edge should exist in G";
	EXPECT_FALSE(G.hasEdge(0, 2)) << "edge should not exist in G";

	G.removeEdge(0, 1);
	EXPECT_FALSE(G.hasEdge(0, 1)) << "edge should no longer exist in G";

}


TEST_F(GraphGTest, testGraphCopy) {
	int64_t n = 3;
	Graph G1(n);
	G1.addEdge(0, 1);

	Graph G2 = G1;	// copy G1 to G2

	EXPECT_EQ(G1.numberOfNodes(), G2.numberOfNodes()) << "G2 should have same number of nodes as G1";
	EXPECT_EQ(G1.numberOfEdges(), G2.numberOfEdges()) << "G2 should have same number of edges as G1";

	G1.addEdge(0, 2);

	EXPECT_FALSE(G2.hasEdge(0, 2)) << "edge inserted in G1 should not exist in G2";
}

TEST_F(GraphGTest, testSubgraphPartitioning) {
	GraphGenerator graphGenerator;
	Graph G = graphGenerator.makeCompleteGraph(4);
	node u = 0;
	node v = 1;
	node w = 2;
	std::unordered_set<node> subgraphSet = {u,v,w};

	Graph subG = Subgraph::fromNodes(G,subgraphSet);
	EXPECT_EQ(3u, subG.numberOfEdges());
	EXPECT_EQ(3u, subG.numberOfNodes());

}


TEST_F(GraphGTest, testBFSIterator) {
	Graph G(6);
	G.addEdge(0,1);
	G.addEdge(1,4);
	G.addEdge(4,5);
	G.addEdge(0,2);
	G.addEdge(0,3);

	std::vector<node> visited;
	G.BFSfrom(0, [&](node v){
		visited.push_back(v);
	});

	std::vector<node> expected = {0,1,2,3,4,5};
	EXPECT_EQ(expected, visited);
}


TEST_F(GraphGTest, testDFSIterator) {
	Graph G(7);
	G.addEdge(0,1);
	G.addEdge(1,4);
	G.addEdge(4,5);
	G.addEdge(0,2);
	G.addEdge(0,3);
	G.addEdge(3,6);

	std::vector<node> visited;
	G.DFSfrom(0, [&](node v){
		visited.push_back(v);
	});

	std::vector<node> expected = {0,3,6,2,1,4,5};
	EXPECT_EQ(expected, visited);
}

TEST_F(GraphGTest, testBFSClass) {
	Graph G(5);
	G.addEdge(0,1);
	G.addEdge(1,2);
	G.addEdge(2,3);
	G.addEdge(3,4);
	BFS bfs(G,1);
	bfs.run();
	auto dist = bfs.getDistances();
	std::vector<count> correct_dist = {1,0,1,2,3};
	for( index i = 0, end = dist.size(); i < end; ++i ) {
		EXPECT_EQ(correct_dist[i],dist[i]) << "distance for node "<<i<<" is not correct";
	}
	auto path = bfs.getPath(4);
	std::vector<node> correct_path = {1,2,3,4};
	for( index i = 0, end = path.size(); i < end; ++i ) {
		TRACE(path);
		EXPECT_EQ(correct_path[i],path[i]) << "node "<<i<<" in path from 1 to 4 is not correct";
	}
}

TEST_F(GraphGTest, testDijkstraClass) {
	Graph G(5,true);
	G.addEdge(0,1,2);
	G.addEdge(1,2,3);
	G.addEdge(2,3,4);
	G.addEdge(3,4,5);
	Dijkstra dij(G,1);
	dij.run();
	auto dist = dij.getDistances();
	std::vector<count> correct_dist = {2,0,3,7,12};
	for( index i = 0, end = dist.size(); i < end; ++i ) {
		EXPECT_EQ(correct_dist[i],dist[i]) << "distance for node "<<i<<" is not correct";
	}
	auto path = dij.getPath(4);
	std::vector<node> correct_path = {1,2,3,4};
	for( index i = 0, end = path.size(); i < end; ++i ) {
		EXPECT_EQ(correct_path[i],path[i]) << "node "<<i<<" in path from 1 to 4 is not correct";
	}
	G.addEdge(0,4,7);
	Dijkstra dij2(G,1);
	dij2.run();
	dist = dij2.getDistances();
	correct_dist = {2,0,3,7,9};
	for( index i = 0, end = dist.size(); i < end; ++i ) {
		EXPECT_EQ(correct_dist[i],dist[i]) << "distance for node "<<i<<" is not correct";
	}
	path = dij2.getPath(4);
	correct_path = {1,0,4};
	for( index i = 0, end = path.size(); i < end; ++i ) {
		EXPECT_EQ(correct_path[i],path[i]) << "node "<<i<<" in path from 1 to 4 is not correct";
	}

}

//TEST_F(GraphGTest, testParallelEdgeInsertion) {
//	count n = 100;
//	Graph G(n);
//
//	#pragma omp parallel for
//	for (node u = 0; u < n; u++) {
//		for (node v = 0; v < n; v++) {
//			if (u < v) {
//				G.addEdge(u, v);
//			}
//		}
//	}
//
//	count expm = ((n - 1) * n) / 2;
//	count m = G.numberOfEdges();
//	EXPECT_EQ(expm, m);
//
//	count d = G.degree(0);
//	EXPECT_EQ((n - 1), d);
//
//}

} /* namespace NetworKit */


void NetworKit::GraphGTest::SetUp() {
}

void NetworKit::GraphGTest::TearDown() {
}

#endif /*NOGTEST */
