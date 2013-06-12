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
	Graph G(0); // empty graphw
	// 2. create proxy
	GraphEventProxy Gproxy(G);
	// 3. create generator and pass proxy
	DynamicGraphGenerator* gen = new DynamicBarabasiAlbertGenerator(Gproxy, 2);
	// 4. create dynamic algorithm and pass graph
	DynamicCommunityDetector* dynPLP = new DynamicLabelPropagation(G, 0, "reactivate");
	// 5. register dynamic algorithm as observer
	Gproxy.registerObserver(dynPLP);
	// 6. initialize graph for generator
	gen->initializeGraph();


	// 7. start generator
	count n1 = 1000;
	count n2 = 2000;

	gen->generateNodes(n1); // stops when graph has n1 nodes
	// 8. start clusterer
	Clustering zeta1 = dynPLP->run();

	EXPECT_TRUE(zeta1.isProper(G)) << "first dynamic clustering should be a proper clustering of G";

	// 9. resume generator
	gen->generateNodes(n2); // terminate when function returns true
	EXPECT_EQ(n2, G.numberOfNodes()) << n2 << "nodes should have been generated";
	// 10. resume clusterer
	Clustering zeta2 = dynPLP->run();
	EXPECT_TRUE(zeta2.isProper(G)) << "second dynamic clustering should be a proper clustering of G";


	INFO("number of clusters 1: " << zeta1.numberOfClusters());
	INFO("number of clusters 2: " << zeta2.numberOfClusters());


	LabelPropagation PLP;
	Clustering zetaPLP = PLP.run(G);
	INFO("number of clusters for static PLP: " << zetaPLP.numberOfClusters());
	EXPECT_TRUE(zetaPLP.isProper(G));

	Louvain PLM;
	Clustering zetaPLM = PLM.run(G);
	INFO("number of clusters for static PLM: " << zetaPLM.numberOfClusters());
	EXPECT_TRUE(zetaPLM.isProper(G));



}

} /* namespace NetworKit */
