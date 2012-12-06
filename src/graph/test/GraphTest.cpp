/*
 * GraphTest.cpp
 *
 *  Created on: 06.12.2012
 *      Author: cls
 */

#include "GraphTest.h"

namespace EnsembleClustering {

CPPUNIT_TEST_SUITE_REGISTRATION(GraphTest);

GraphTest::GraphTest() {
	// TODO Auto-generated constructor stub

}

GraphTest::~GraphTest() {
	// TODO Auto-generated destructor stub
}

void GraphTest::setUp() {
	this->randomGraph = this->gen.makeErdosRenyiGraph(20, 0.2);
}

void GraphTest::tearDown() {
}


void GraphTest::testIteration() {
	INFO("testing iteration");
	int success = 0;

	Graph G = this->randomGraph;

	int64_t etype = G.defaultEdgeType;

	STINGER_PARALLEL_FORALL_EDGES_BEGIN(G.asSTINGER(), etype) {
		node u = STINGER_EDGE_SOURCE;
		node v = STINGER_EDGE_DEST;
		TRACE("found edge (" << u << "," << v << ") with weight " << stinger_edgeweight(G.asSTINGER(), u, v, etype));
	} STINGER_PARALLEL_FORALL_EDGES_END();

	success = 1;
	CPPUNIT_ASSERT_EQUAL(success, 1);

}

void GraphTest::testCppUnit() {
	int x = 1;
	CPPUNIT_ASSERT_EQUAL(x, 1);

}

} /* namespace EnsembleClustering */
