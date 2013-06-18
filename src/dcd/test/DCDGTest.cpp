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
	DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(2);

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

	DynamicGraphSource* dynGen = new DynamicPubWebGenerator(numInitialNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors);

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
	DynamicGraphSource* dynGen = new DynamicBarabasiAlbertGenerator(k);

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
	 DynamicGraphSource* dynGen = new DynamicBarabasiAlbertGenerator(1);
	 DynamicCommunityDetector* dynCD1 = new DynamicLabelPropagation(0, "Reactivate");
	 DynamicCommunityDetector* dynCD2 = new DynamicLabelPropagation(0, "ReactivateNeighbors");

	 std::vector<DynamicCommunityDetector*> detectors = {dynCD1, dynCD2};
	 DynCDSetup setup(*dynGen, detectors, 10e4, 100);

	 setup.run();

}


TEST_F(DCDGTest, tryDGSAsSource) {
	std::string path;
	std::cout << "[INPUT] .dgs file path >" << std::endl;
	std::getline(std::cin, path);

	DynamicGraphSource* source = new DynamicDGSParser(path);
	DynamicCommunityDetector* dynCD1 = new DynamicLabelPropagation(0, "Reactivate");
	// DynamicCommunityDetector* dynCD2 = new DynamicLabelPropagation(0, "ReactivateNeighbors");

	std::vector<DynamicCommunityDetector*> detectors = { dynCD1 };
	DynCDSetup setup(*source, detectors, 10e9, 10000);

	setup.run();
}


TEST_F(DCDGTest, testPseudoDynamic) {
	Graph G(3);
	G.addEdge(0, 1);
	G.addEdge(1, 2);


	DynamicGraphSource* source = new PseudoDynamic(G);

	GraphEventProxy* proxy = source->newGraph();
	Graph* dynG = proxy->G;

	source->initializeGraph();

	source->generate();
	EXPECT_EQ(1, dynG->numberOfNodes());
	EXPECT_EQ(0, dynG->numberOfEdges());

	source->generate();
	EXPECT_EQ(2, dynG->numberOfNodes());
	EXPECT_EQ(1, dynG->numberOfEdges());


	source->generate();
	EXPECT_EQ(3, dynG->numberOfNodes());
	EXPECT_EQ(2, dynG->numberOfEdges());


	EXPECT_EQ(G.numberOfNodes(), dynG->numberOfNodes()) << "n should be equal";
	EXPECT_EQ(G.numberOfEdges(), dynG->numberOfEdges()) << "m should be equal";

}

TEST_F(DCDGTest, testPseudoDynamicOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");

	DynamicGraphSource* source = new PseudoDynamic(G);

	Graph* dynG = source->newGraph()->G;

	source->initializeGraph();

	source->generateNodes(G.numberOfNodes());

	EXPECT_EQ(G.numberOfEdges(), dynG->numberOfEdges()) << "same graph, so m should be equal";
}


} /* namespace NetworKit */

#endif /*NOGTEST */
