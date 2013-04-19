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


TEST_F(DynamicCommunityDetectionGTest, testDynamicLabelPropagation) {
	// 1. create graph
	Graph G(0); // empty graph
	// 2. create proxy
	GraphEventProxy Gproxy(G);
	// 3. create generator and pass proxy
	DynamicGraphGenerator* gen = new DynamicBarabasiAlbertGenerator(Gproxy, 3);
	// 4. create dynamic algorithm and pass graph
	DynamicClusterer* dynPLP = new DynamicLabelPropagation(G, 0, "reactivate");
	// 5. register dynamic algorithm as observer
	Gproxy.registerObserver(dynPLP);
	// 6. initialize graph for generator
	gen->initializeGraph();


	// 7. start generator
	count n1 = 1000;
	count n2 = 2000;

	auto cont1 = [&](){
		return (G.numberOfNodes() < n1);
	};
	gen->generateWhile(cont1); // stops when function returns true
	// 8. start clusterer
	Clustering zeta1 = dynPLP->run();
	// 9. resume generator
	auto cont2 = [&](){
		return (G.numberOfNodes() < n2);
	};
	gen->generateWhile(cont2); // terminate when function returns true
	EXPECT_EQ(n2, G.numberOfNodes()) << n2 << "nodes should have been generated";
	// 10. resume clusterer
	Clustering zeta2 = dynPLP->run();

	INFO("number of clusters 1: " << zeta1.numberOfClusters());
	INFO("number of clusters 2: " << zeta2.numberOfClusters());


	LabelPropagation PLP;
	Clustering zetaStatic = PLP.run(G);
	INFO("number of clusters static: " << zetaStatic.numberOfClusters());


}

} /* namespace NetworKit */
