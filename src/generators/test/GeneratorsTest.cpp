/*
 * GeneratorsTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "GeneratorsTest.h"

namespace NetworKit {

GeneratorsTest::GeneratorsTest() {
	// TODO Auto-generated constructor stub

}

GeneratorsTest::~GeneratorsTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(GeneratorsTest, testDynamicBarabasiAlbertGenerator) {

	Graph G(0); // empty graph
	GraphEventProxy proxy(G);

	DynamicGraphGenerator* gen = new DynamicBarabasiAlbertGenerator(proxy);

	gen->initializeGraph();

	EXPECT_EQ(2, G.numberOfNodes()) << "initially the generator creates two connected nodes";

	count n = 100;

	gen->generate([&](){
		return ( G.numberOfNodes() == n );
	});

	EXPECT_EQ(n, G.numberOfNodes());

	DEBUG("m = " << G.numberOfEdges());
}

} /* namespace NetworKit */
