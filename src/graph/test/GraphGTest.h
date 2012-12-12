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

} /* namespace EnsembleClustering */
#endif /* GRAPHGTEST_H_ */
