/*
 * GraphGTest.h
 *
 *  Created on: 12.12.2012
 *      Author: cls
 */

#ifndef GRAPHGTEST_H_
#define GRAPHGTEST_H_

#include <gtest/gtest.h>

#include "../../aux/log.h"
#include "../Graph.h"
#include "../GraphGenerator.h"


namespace EnsembleClustering {

class GraphGTest: public testing::Test {

protected:

	GraphGenerator gen;

public:

	virtual void SetUp();

	virtual void TearDown();

};

TEST_F(GraphGTest, testIteration) {
	int success = 0;

	Graph G = this->gen.makeCompleteGraph(20);

	int64_t etype = G.defaultEdgeType;

	STINGER_PARALLEL_FORALL_EDGES_BEGIN(G.asSTINGER(), etype) {
		node u = STINGER_EDGE_SOURCE;
		node v = STINGER_EDGE_DEST;
		TRACE("found edge (" << u << "," << v << ") with weight " << stinger_edgeweight(G.asSTINGER(), u, v, etype));
	} STINGER_PARALLEL_FORALL_EDGES_END();

	success = 1;
	EXPECT_EQ(success, 1);
}


TEST_F(GraphGTest, testUndirectedEdgeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t edgeCount = 0;
	READ_ONLY_FORALL_EDGES_BEGIN(G) {
		node u = EDGE_SOURCE;
		node v = EDGE_DEST;
		if (u < v) {
			edgeCount += 1;
		}
	} READ_ONLY_FORALL_EDGES_END();

	EXPECT_EQ(edgeCount, (n * (n-1)) / 2) << "There are (n * (n-1)) / 2 undirected edges in a compete graph";
}


TEST_F(GraphGTest, testLambdaEdgeIteration) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	int64_t edgeCount = 0;
	G.forallEdges([&](node u, node v) {
		if (u < v) {
			edgeCount += 1;
		}
	});

	EXPECT_EQ(edgeCount, (n * (n-1)) / 2) << "There are (n * (n-1)) / 2 undirected edges in a compete graph";
}


TEST_F(GraphGTest, testLambdaEdgeModification) {

	int64_t n = 100;
	Graph G = this->gen.makeCompleteGraph(n);

	G.forallEdges([&](node u, node v){
		G.removeEdge(u, v);
	});

	EXPECT_EQ(G.numberOfEdges(), 0) << "all edges should have been deleted";
}


} /* namespace EnsembleClustering */
#endif /* GRAPHGTEST_H_ */
