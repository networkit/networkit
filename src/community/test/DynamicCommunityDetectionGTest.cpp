/*
 * DynamicCommunityDetectionGTest.cpp
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#include "DynamicCommunityDetectionGTest.h"

namespace NetworKit {

DynamicCommunityDetectionGTest::DynamicCommunityDetectionGTest() {
	// TODO Auto-generated constructor stub

}

DynamicCommunityDetectionGTest::~DynamicCommunityDetectionGTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(DynamicCommunityDetectionGTest, testSetup) {
	// 1. create graph
	Graph G(0); // empty graph
	// 2. create proxy
	GraphEventProxy Gproxy(G);
	// 3. create generator and pass proxy
	DynamicGraphGenerator* gen = new DynamicBarabasiAlbertGenerator(Gproxy);
	// 4. create dynamic algorithm and pass graph
	DynamicClusterer* dynPLP = new DynamicLabelPropagation(G);
	// 5. register dynamic algorithm as observer
	Gproxy.registerObserver(dynPLP);
	// 6. start generator
	auto until1 = [&](){
		return (G.numberOfNodes() == 100);
	};
	gen->generate(until1); // stops when function returns true
	// 7. start clusterer
	Clustering zeta1 = dynPLP->run();
	// 8. resume generator
	auto until2 = [&](){
		return (G.numberOfNodes() == 200);
	};
	gen->generate(until2); // terminate when function returns true
	EXPECT_EQ(200, G.numberOfNodes()) << "200 nodes should have been generated";
	// 9. resume clusterer
	Clustering zeta2 = dynPLP->run();

}

} /* namespace NetworKit */
