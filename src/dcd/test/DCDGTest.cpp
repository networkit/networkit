/*
 * DynamicCommunityDetectionGTest.cpp
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#ifndef NOGTEST


#include "DCDGTest.h"

namespace NetworKit {

DCDGTest::DCDGTest() {
	// TODO Auto-generated constructor stub

}

DCDGTest::~DCDGTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(DCDGTest, testDynamicLabelPropagation) {

	//  create generator and pass proxy
	DynamicGraphGenerator* gen = new DynamicBarabasiAlbertGenerator(2);

	GraphEventProxy* Gproxy = gen->newGraph();
	Graph* G = Gproxy->G;

	// create dynamic algorithm and pass graph
	DynamicCommunityDetector* dynCD = new DynamicLabelPropagation(0, "Reactivate");
	dynCD->setGraph(*G);
	// 5. register dynamic algorithm as observer
	Gproxy->registerObserver(dynCD);
	// 6. initialize graph for generator
	gen->initializeGraph();


	// 7. start generator
	count n1 = 1000;
	count n2 = 2000;

	gen->generateNodes(n1); // stops when graph has n1 nodes
	// 8. start clusterer
	Clustering zeta1 = dynCD->run();

	EXPECT_TRUE(zeta1.isProper(*G)) << "first dynamic clustering should be a proper clustering of G";

	// 9. resume generator
	gen->generateNodes(n2); // terminate when function returns true
	EXPECT_EQ(n2, G->numberOfNodes()) << n2 << "nodes should have been generated";
	// 10. resume clusterer
	Clustering zeta2 = dynCD->run();
	EXPECT_TRUE(zeta2.isProper(*G)) << "second dynamic clustering should be a proper clustering of G";


	INFO("number of clusters 1: " << zeta1.numberOfClusters());
	INFO("number of clusters 2: " << zeta2.numberOfClusters());


	LabelPropagation PLP;
	Clustering zetaPLP = PLP.run(*G);
	INFO("number of clusters for static PLP: " << zetaPLP.numberOfClusters());
	EXPECT_TRUE(zetaPLP.isProper(*G));

	Louvain PLM;
	Clustering zetaPLM = PLM.run(*G);
	INFO("number of clusters for static PLM: " << zetaPLM.numberOfClusters());
	EXPECT_TRUE(zetaPLM.isProper(*G));

}


TEST_F(DCDGTest, tryArxivGraphs) {
	std::string graphPath;
	std::cout << "[INPUT] .dgs file path >" << std::endl;
	std::getline(std::cin, graphPath);


}


TEST_F(DCDGTest, tryDynamicPubWebGeneratorAsSource) {
	count numInitialNodes = 600;
	count numberOfDenseAreas = 10;
	float neighborhoodRadius = 0.125;
	count maxNumberOfNeighbors = 16;
	count numIterations = 10;

	DynamicGraphGenerator* dynGen = new DynamicPubWebGenerator(numInitialNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors);

	GraphEventProxy* Gproxy = dynGen->newGraph();
	Graph* G = Gproxy->G;

	DynamicCommunityDetector* dynCD = new DynamicLabelPropagation(0, "Reactivate");
	dynCD->setGraph(*G);




	Gproxy->registerObserver(dynCD);

	count deltaT = 1;
	count tMax = 10;

	dynGen->initializeGraph();

	std::vector<Clustering> results;
	while (G->time() <= tMax) {
		dynGen->generateTimeSteps(G->time() + deltaT);
		if (G->time() % 2 == 0) {
			results.push_back(dynCD->run());
		}
	}

	for (Clustering zeta : results) {
		DEBUG("number of clusters: " << zeta.numberOfClusters());
	}

}


TEST_F(DCDGTest, tryDynamicBarabasiAlbertGeneratorAsSource) {

	// create instance of generator
	count k = 1;
	DynamicGraphGenerator* dynGen = new DynamicBarabasiAlbertGenerator(k);

	// use generator to get graph and proxy
	GraphEventProxy* Gproxy = dynGen->newGraph();
	Graph* G = Gproxy->G;

	// create instance of dynamic community detection algorithm
	DynamicCommunityDetector* dynCD = new DynamicLabelPropagation(0, "ReactivateNeighbors");

	// provide graph to algorithm
	dynCD->setGraph(*G);
	// register algorithm as observer to the graph
	Gproxy->registerObserver(dynCD);

	// initialize graph
	dynGen->initializeGraph();

	// start dynamic community detection setup
	count deltaT = 1;
	count tMax = 10;

	std::vector<Clustering> results;
	while (G->time() <= tMax) {
		dynGen->generateTimeSteps(G->time() + deltaT);
		if (G->time() % 2 == 0) {
			results.push_back(dynCD->run());
		}
	}

	for (Clustering zeta : results) {
		DEBUG("number of clusters: " << zeta.numberOfClusters());
	}

}


TEST_F(DCDGTest, tryDynCDSetup) {
	 DynamicGraphGenerator* dynGen = new DynamicBarabasiAlbertGenerator(1);
	 DynamicCommunityDetector* dynCD1 = new DynamicLabelPropagation(0, "Reactivate");
	 DynamicCommunityDetector* dynCD2 = new DynamicLabelPropagation(0, "ReactivateNeighbors");

	 std::vector<DynamicCommunityDetector*> detectors = {dynCD1, dynCD2};
	 DynCDSetup setup(*dynGen, detectors, 10e4, 100);

	 setup.run();

}

} /* namespace NetworKit */

#endif /*NOGTEST */
