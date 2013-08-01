/*
 * Graph2GTest.cpp
 *
 *  Created on: 04.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "Graph2GTest.h"

namespace NetworKit {

Graph2GTest::Graph2GTest() {
	// TODO Auto-generated constructor stub

}

Graph2GTest::~Graph2GTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(Graph2GTest, testNumberOfNodes) {
	int64_t n = 10;
	Graph G(n);

	EXPECT_EQ(n, G.numberOfNodes()) << "number of nodes should be " << n << " in the graph";


	// TODO: insert nodes and delete them, test again
}

TEST_F(Graph2GTest, testNumberOfEdges) {
	int64_t n = 10;
	Graph G(n);

	EXPECT_EQ(0, G.numberOfEdges()) << "no edges should be in the graph";

	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(1, 3);

	EXPECT_TRUE(G.hasEdge(1, 2));
	EXPECT_TRUE(G.hasEdge(2, 3));
	EXPECT_TRUE(G.hasEdge(1, 3));
	EXPECT_EQ(3, G.numberOfEdges()) << "3 edges should have been inserted";

	G.removeEdge(3, 1);
	G.removeEdge(2, 1);
	G.removeEdge(3, 2);

	EXPECT_FALSE(G.hasEdge(1, 2));
	EXPECT_FALSE(G.hasEdge(2, 3));
	EXPECT_FALSE(G.hasEdge(1, 3));
	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

TEST_F(Graph2GTest, testIsEmpty) {
	Graph G(0);
	EXPECT_EQ(0, G.numberOfNodes()) << "no nodes should be in the graph";
	EXPECT_TRUE(G.isEmpty());

	// TODO: insert node, check if isEmpty is false, delete nodes again, check if true again
}

TEST_F(Graph2GTest, testEdgeInsertionAndRemoval) {

	int64_t n = 10;
	Graph G(n);

	EXPECT_EQ(0, G.numberOfEdges()) << "no edges should be in the graph";

	G.addEdge(1, 2);
	G.addEdge(2, 3);
	G.addEdge(1, 3);

	EXPECT_TRUE(G.hasEdge(1, 2));
	EXPECT_TRUE(G.hasEdge(2, 3));
	EXPECT_TRUE(G.hasEdge(1, 3));
	EXPECT_EQ(3, G.numberOfEdges()) << "3 edges should have been inserted";

	G.removeEdge(3, 1);
	G.removeEdge(2, 1);
	G.removeEdge(3, 2);

	EXPECT_FALSE(G.hasEdge(1, 2));
	EXPECT_FALSE(G.hasEdge(2, 3));
	EXPECT_FALSE(G.hasEdge(1, 3));
	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

TEST_F(Graph2GTest, testNodeIteration) {
	int64_t n = 20;
	Graph G(n);

	G.forNodes([&](node v) {
		G.addEdge(v, (v+1) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	G.forNodes([&](node v) {
		G.removeEdge(v, (v+1) % n);
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}


//TEST_F(Graph2GTest, testParallelEdgeInsertion) {
//	EXPECT_TRUE(false) << "TODO: implement";
//	return;
//
//
//	int64_t n = 500;
//	int64_t offset = 100;
//	Graph G(n);
//
//	try {
//		G.parallelForNodes([&](node v) {
//			G.addEdge(v, (v+offset) % n);
//		});
//
//	} catch (...) {
//		EXPECT_TRUE(false) << "exception thrown - fail test";
//	}
//
//	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";
//
//	G.forNodes([&](node v) {
//		EXPECT_EQ(2, G.degree(v)) << "degree should be two";
//	});
//
//}


//TEST_F(Graph2GTest, testParallelEdgeRemoval) {
//	EXPECT_TRUE(false) << "TODO: implement";
//	return;
//
//	int64_t n = 500;
//	int64_t offset = 100;
//	Graph G(n);
//
//	G.forNodes([&](node v) {
//		G.addEdge(v, (v+offset) % n);
//	});
//
//	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";
//
//	try {
//		G.parallelForNodes([&](node v) {
//			G.removeEdge(v, (v+offset) % n);
//		});
//	} catch (...) {
//		EXPECT_TRUE(false) << "exception thrown - fail test";
//	}
//
//
//	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
//}

TEST_F(Graph2GTest, testEdgeIteration) {
	int64_t n = 100;
	int64_t offset = 20;
	Graph G(n);

	G.forNodes([&](node v) {
		G.addEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	G.forEdges([&](node u, node v) {
		DEBUG("about to remove edge " << u << ", " << v);
		G.removeEdge(u, v);
		DEBUG("have removed edge " << u << ", " << v);
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

//TEST_F(Graph2GTest, testParallelEdgeModification) {
//	EXPECT_TRUE(false) << "TODO: implement";
//	return;
//
//
//	int64_t n = 500;
//	int64_t offset = 100;
//	Graph G(n);
//
//	G.parallelForNodes([&](node v) {
//		G.addEdge(v, (v+offset) % n);
//	});
//
//	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";
//
//	G.parallelForEdges([&](node u, node v) {
//		G.removeEdge(u, v);
//	});
//
//	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
//}

TEST_F(Graph2GTest, testNeighborIteration) {
	int64_t n = 50;
	int64_t offset = 10;
	Graph G(n);

	G.forNodes([&](node v) {
		G.addEdge(v, (v+offset) % n);
		G.addEdge(v, (v+2*offset) % n);
	});

	G.forNodes([&](node v) {
		G.forNeighborsOf(v, [&](node u) {
			EXPECT_TRUE(G.hasEdge(v, u));
		});
	});
}
/*
TEST_F(Graph2GTest, testParallelSumForNodes) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.parallelForNodes([&](node v) {
		G.addEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	// sum degrees in parallel
	float sum;
	float summe = G.parallelSumForNodes([&](node v) {
		return (float) G.degree(v);
	});

	EXPECT_EQ(summe, 2 * G.numberOfEdges()) << " graph volume should be" << (2*G.numberOfEdges());
}
*/
TEST_F(Graph2GTest, testNodePairIteration) {
	int64_t n = 50;
	Graph G(n);

	G.forNodePairs([&](node u, node v) {
		G.addEdge(u, v);
	});

	EXPECT_EQ(n*(n-1)/2, G.numberOfEdges()) << n*(n-1)/2 << " edges should have been inserted";

	G.forNodePairs([&](node u, node v) {
		G.removeEdge(u, v);
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

TEST_F(Graph2GTest, testConstNodeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G2(n);

	G2.forNodes([&](node v) {
		G2.addEdge(v, (v+offset) % n);
	});

	const Graph G(G2);
	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";
	G.forNodes([&](node v) {
		EXPECT_EQ(2, G.degree(v)) << "degree should be two";
	});
}

/*
TEST_F(Graph2GTest, testConstParallelNodeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G2(n);

	G2.parallelForNodes([&](node v) {
		G2.addEdge(v, (v+offset) % n);
	});

	const Graph G(G2);
	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";
	G.parallelForNodes([&](node v) {
		EXPECT_EQ(2, G.degree(v)) << "degree should be two";
	});
}
*/


TEST_F(Graph2GTest, testConstEdgeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.forNodes([&](node v) {
		G.addEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	const Graph G2(G);
	G2.forEdges([&](node u, node v) {
		EXPECT_TRUE(G2.hasEdge(v, u));
	});
}

// Will fail, since the data structure is not thread-safe
/*
TEST_F(Graph2GTest, testConstParallelEdgeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.parallelForNodes([&](node v) {
		G.addEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	const Graph G2(G);
	G2.parallelForEdges([&](node u, node v) {
		EXPECT_TRUE(G2.hasEdge(v, u));
	});
}
*/

TEST_F(Graph2GTest, testConstNeighborIteration) {
	int64_t n = 50;
	int64_t offset = 10;
	Graph G(n);

	G.forNodes([&](node v) {
		G.addEdge(v, (v+offset) % n);
		G.addEdge(v, (v+2*offset) % n);
	});

	const Graph G2(G);
	G2.forNodes([&](node v) {
		G2.forNeighborsOf(v, [&](node u) {
			EXPECT_TRUE(G.hasEdge(v, u));
		});
	});
}
/*
TEST_F(Graph2GTest, testConstParallelSumForNodes) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.parallelForNodes([&](node v) {
		node w = (v + offset) % n;
		assert (w < n);
		G.addEdge(v, w);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	const Graph G2(G);
	// sum degrees in parallel
	float sum;
	float summe = G2.parallelSumForNodes([&](node v) {
		return (float) G.degree(v);
	});

	EXPECT_EQ(summe, 2 * G2.numberOfEdges()) << " graph volume should be" << (2*G2.numberOfEdges());
}*/

TEST_F(Graph2GTest, testConstNodePairIteration) {
	int64_t n = 50;
	Graph G(n);

	G.forNodePairs([&](node u, node v) {
		G.addEdge(u, v);
	});

	EXPECT_EQ(n*(n-1)/2, G.numberOfEdges()) << n*(n-1)/2 << " edges should have been inserted";

	const Graph G2(G);
	G2.forNodePairs([&](node u, node v) {
		EXPECT_TRUE(G2.hasEdge(u, v));
	});
}


TEST_F(Graph2GTest, testAddNode) {
	count n = 10;
	count newN = n + 1;
	Graph G(n);

	node v = G.addNode();
	EXPECT_EQ(v, n) << "new node should have id equal to previous n";
	EXPECT_EQ(newN, G.numberOfNodes()) << "new number of nodes should be " << newN;
	EXPECT_EQ(0, G.degree(v)) << "new node should have degree 0";

}


TEST_F(Graph2GTest, testRemoveNode) {

	// test removing all nodes sequentially
	count n = 10;
	Graph G(n);

	EXPECT_EQ(n, G.numberOfNodes());

	count nc = 0;
	G.forNodes([&](node u){
		nc += 1;
	});
	EXPECT_EQ(nc, n);

	G.forNodes([&](node u) {
		G.removeNode(u);
	});

	EXPECT_EQ(0, G.numberOfNodes()) << "number of nodes should be zero";

	nc = 0;
	G.forNodes([&](node u){
		nc += 1;
	});
	EXPECT_EQ(0, nc) << "counter should not be incremented";


	// test removing a single node
	Graph G2(n);
	node x = 5;

	G.forNodes([&](node u){
		G2.addEdge(x, u);
	});

	G2.removeNode(x);

	G2.forNodes([&](node u) {
		EXPECT_NE(x, u) << "node x should not be handled";
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed with x";


}


TEST_F(Graph2GTest, testWeightedDegree) {
	count n = 4;
	Graph G(n);

	node v = 0;
	G.addEdge(v, 1, 20.0);
	G.addEdge(v, 2, 20.0);
	G.addEdge(v, v, 2.0); // self-loop
	G.addEdge(v, 3, 100.0);
	G.removeEdge(v, 3);

	edgeweight wDeg = G.weightedDegree(v);

	EXPECT_EQ(42.0, wDeg) << "weighted degree should be 42.0";

}


TEST_F(Graph2GTest, testNodeDynamics) {
	count n = 10;
	Graph G(n);

	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});

	for (count i = 0; i < n; ++i) {
		G.addNode();
	}

	EXPECT_EQ((n * 2), G.numberOfNodes());
	EXPECT_EQ((n * (n-1)) / 2, G.numberOfEdges());

	for (node u = 0; u < n; ++u) {
		G.forEdgesOf(u, [&](node u, node v){
			G.removeEdge(u, v);
		});
		G.removeNode(u);
	}

	EXPECT_EQ(n, G.numberOfNodes());
	EXPECT_EQ(0, G.numberOfEdges());

	for (node u = n; u < 2 * n; ++u) {
		G.forEdgesOf(u, [&](node u, node v){
			G.removeEdge(u, v);
		});
		G.removeNode(u);
	}

	EXPECT_EQ(0, G.numberOfNodes());

}


TEST_F(Graph2GTest, testSetWeight) {
	count n = 4;
	Graph G(n);

	node u = 0;
	node v = 1;
	G.addEdge(u, v);

	G.setWeight(u, v, 42.0);
	EXPECT_EQ(42.0, G.weight(u, v));

	// symmetric case
	G.setWeight(v, u, 23.0);
	EXPECT_EQ(23.0, G.weight(u, v));

	// self-loop
	G.addEdge(v, v);
	G.setWeight(v, v, 42.0);
	EXPECT_EQ(42.0, G.weight(v, v));
	G.setWeight(v, v, 17.0);
	EXPECT_EQ(17.0, G.weight(v, v));


}


TEST_F(Graph2GTest, testForWeightedEdges) {
	count n = 4;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v, 1.0);
	});

	edgeweight weightSum = 0.0;
	G.forWeightedEdges([&](node u, node v, node w){
		weightSum += w;
	});

	EXPECT_EQ(6.0, weightSum) << "sum of edge weights should be 6";

}


TEST_F(Graph2GTest, testParallelForWeightedEdges) {
	count n = 4;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v, 1.0);
	});

	edgeweight weightSum = 0.0;
	G.parallelForWeightedEdges([&](node u, node v, node w){
		#pragma omp atomic update
		weightSum += w;
	});

	EXPECT_EQ(6.0, weightSum) << "sum of edge weights should be 6";

}


TEST_F(Graph2GTest, testEdgeAttributes) {
	count n = 5;
	Graph G(n);

	int attrId = G.addEdgeAttribute_double(0.0);

	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});

	G.forEdges([&](node u, node v){
		EXPECT_EQ(0.0, G.attribute_double(u, v, attrId));
	});

	G.forEdges([&](node u, node v){
		G.setAttribute_double(u, v, attrId, 42.0);
	});

	G.forEdges([&](node u, node v){
		EXPECT_EQ(42.0, G.attribute_double(u, v, attrId));
	});

	node v = G.addNode();
	G.addEdge(v, 0);

	EXPECT_EQ(0.0, G.attribute_double(v, 0, attrId));

}


TEST_F(Graph2GTest, testTotalEdgeWeight) {
	count n = 100;
	Graph G(n);

	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});

	double m = (n * (n-1)) / 2.0;
	double total = G.totalEdgeWeight();
	EXPECT_EQ(m, total) << "total edge weight should be n * (n -1) / 2";

	INFO("difference: " << (total - m));

}


TEST_F(Graph2GTest, testNodeAdditionAndIteration) {
	count n = 1;
	Graph G(n);

	G.addNode();

	EXPECT_EQ((n + 1), G.numberOfNodes());

	count nc = 0;
	G.forNodes([&](node u){
		nc += 1;
	});

	EXPECT_EQ((n + 1), nc);
}


TEST_F(Graph2GTest, testConditionalNodeIteration) {
	count n = 100;
	Graph G(n);

	count nc = 0;
	G.forNodesWhile([&](){ return (nc < n); }, [&](node u){
		nc += 1;
	});

	EXPECT_EQ(n, nc);
}


TEST_F(Graph2GTest, testMultiEdgeInsertion) {
	// the current implementation of Graph allows multi-edges for performance reasons
	// this may result in unexpected behavior - it is up to the graph generator to avoid multi-edges if necessary
	count n = 2;
	Graph G(n);

	node u = 0;
	node v = 1;
	for (int i = 0; i < 42; ++i) {
		G.addEdge(u, v);
	}

	EXPECT_EQ(42, G.numberOfEdges());
	EXPECT_EQ(42, G.degree(u));
	EXPECT_EQ(42, G.degree(v));
}


TEST_F(Graph2GTest, testNodeCollection) {
	count n = 4;
	Graph G(n);
	G.removeNode(0);
	std::vector<node> nodes = G.nodes();
	EXPECT_EQ(G.numberOfNodes(), nodes.size());
}

TEST_F(Graph2GTest, testEdgeCollection) {
	count n = 4;
	Graph G(n);
	G.addEdge(0, 1);
	std::vector<std::pair<node, node> > edges = G.edges();
	EXPECT_EQ(G.numberOfEdges(), edges.size());
}


TEST_F(Graph2GTest, testHasNode) {
	Graph G(4);
	EXPECT_TRUE(G.hasNode(3));
	EXPECT_TRUE(G.hasNode(2));
	EXPECT_TRUE(G.hasNode(1));
	EXPECT_TRUE(G.hasNode(0));
}


} /* namespace NetworKit */

#endif /*NOGTEST */
