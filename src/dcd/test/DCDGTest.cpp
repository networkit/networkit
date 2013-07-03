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
	DynamicCommunityDetector* dynCD = new DynamicLabelPropagation(0, "Isolate");
	dynCD->setGraph(*G);
	// 5. register dynamic algorithm as observer
	Gproxy->registerObserver(dynCD);
	// 6. initialize graph for generator
	gen->initializeGraph();


	// 7. start generator
	count n1 = 10;
	count n2 = 20;

	gen->generateNodes(n1); // stops when graph has n1 nodes
	// 8. start clusterer
	Clustering zeta1 = dynCD->run();

	EXPECT_EQ(G->numberOfNodes(), zeta1.numberOfEntries()) << "clustering must have as many entries as nodes";

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

	INFO("first clustering: " << Aux::vectorToString(zeta1.getVector()));

	Louvain PLM;
	Clustering zetaPLM = PLM.run(*G);
	INFO("number of clusters for static PLM: " << zetaPLM.numberOfClusters());
	EXPECT_TRUE(zetaPLM.isProper(*G));

	INFO("second clustering: " << Aux::vectorToString(zeta1.getVector()));

}


TEST_F(DCDGTest, tryArxivGraphs) {
	std::string graphPath;
//	std::cout << "[INPUT] .dgs file path >" << std::endl;
//	std::getline(std::cin, graphPath);

	graphPath = "/Users/cls/workspace/Data/arXiv/CS-all-paper.dgs";

	DynamicGraphSource* source = new DynamicDGSParser(graphPath);

	DynamicCommunityDetector* dynCD = new DynamicLabelPropagation(0, "Isolate");

	std::vector<DynamicCommunityDetector*> targets = { dynCD };
	DynCDSetup setup(*source, targets, 1e9, 1000);

	setup.run();


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
	 DynamicCommunityDetector* dynCD1 = new DynamicLabelPropagation(0, "Isolate");
	 DynamicCommunityDetector* dynCD2 = new DynamicLabelPropagation(0, "IsolateNeighbors");

	 std::vector<DynamicCommunityDetector*> detectors = {dynCD1, dynCD2};
	 DynCDSetup setup(*dynGen, detectors, 1e4, 100);

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


TEST_F(DCDGTest, tryStaticVsDynamic) {
	std::string path;
	std::cout << "[INPUT] .graph file path >" << std::endl;
	std::getline(std::cin, path);

	METISGraphReader reader;

	Graph Gstatic = reader.read(path);
	DynamicGraphSource* source = new PseudoDynamic(Gstatic);
	DynamicCommunityDetector* dynCD1 = new DynamicLabelPropagation(0, "Isolate");

	std::vector<DynamicCommunityDetector*> detectors = { dynCD1 };
	DynCDSetup setup(*source, detectors, 10e9, 10);

	setup.run();

	Graph* G = setup.getGraph();
	LabelPropagation PLP;
	Clustering zetaStatic = PLP.run(*G);

	INFO("number of clusters for static: " << zetaStatic.numberOfClusters());
}


TEST_F(DCDGTest, testPseudoDynamic) {
	Graph G(4);
	G.removeNode(3);
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

	EXPECT_EQ(G.numberOfNodes(), dynG->numberOfNodes()) << "n should be equal";
	EXPECT_EQ(G.numberOfEdges(), dynG->numberOfEdges()) << "m should be equal";
}



TEST_F(DCDGTest, testStrategyReactivate) {
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");
	DynamicGraphSource* dynGen = new PseudoDynamic(G);
	DynamicCommunityDetector* dynCD1 = new DynamicLabelPropagation(0, "Reactivate");

	std::vector<DynamicCommunityDetector*> detectors = { dynCD1 };
	DynCDSetup setup(*dynGen, detectors, 200, 50);

	setup.run();
}


TEST_F(DCDGTest, testStrategyReactivateNeighbors) {
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");
	DynamicGraphSource* dynGen = new PseudoDynamic(G);
	DynamicCommunityDetector* dynCD1 = new DynamicLabelPropagation(0, "ReactivateNeighbors");

	 std::vector<DynamicCommunityDetector*> detectors = {dynCD1};
	 DynCDSetup setup(*dynGen, detectors, 200, 50);

	 setup.run();
}


TEST_F(DCDGTest, testStrategyIsolate) {
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");
	DynamicGraphSource* dynGen = new PseudoDynamic(G);
	DynamicCommunityDetector* dynCD1 = new DynamicLabelPropagation(0, "Isolate");

	 std::vector<DynamicCommunityDetector*> detectors = {dynCD1};
	 DynCDSetup setup(*dynGen, detectors, 200, 50);

	 setup.run();
}


TEST_F(DCDGTest, testStrategyIsolateNeighbors) {
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");
	DynamicGraphSource* dynGen = new PseudoDynamic(G);
	DynamicCommunityDetector* dynCD1 = new DynamicLabelPropagation(0, "IsolateNeighbors");

	 std::vector<DynamicCommunityDetector*> detectors = {dynCD1};
	 DynCDSetup setup(*dynGen, detectors, 200, 50);

	 setup.run();
}


TEST_F(DCDGTest, tryDynamicEnsemble) {
	 DynamicGraphSource* dynGen = new DynamicBarabasiAlbertGenerator(1);
	 DynamicCommunityDetector* dynLP1 = new DynamicLabelPropagation(0, "Isolate");
	 DynamicCommunityDetector* dynLP2 = new DynamicLabelPropagation(0, "IsolateNeighbors");
	 Clusterer* PLM = new Louvain();

	 DynamicEnsemble* ensemble = new DynamicEnsemble();
	 ensemble->addBaseAlgorithm(*dynLP1);
	 ensemble->addBaseAlgorithm(*dynLP2);
	 ensemble->setFinalAlgorithm(*PLM);

	 HashingOverlapper overlapAlgo;
	 ensemble->setOverlapper(overlapAlgo);

	 std::vector<DynamicCommunityDetector*> detectors = {ensemble};
	 DynCDSetup setup(*dynGen, detectors, 1e4, 1000);

	 setup.run();

}

TEST_F(DCDGTest, testDynamicLabelPropagation2) {

	DynamicCommunityDetector* dynPLP = new DynamicLabelPropagation(0, "Isolate");
	INFO("created algorithm: " << dynPLP->toString());

	DynamicGraphSource* dynGen = new DynamicBarabasiAlbertGenerator(1);


	std::vector<DynamicCommunityDetector*> detectors = {dynPLP};
	DynCDSetup setup(*dynGen, detectors, 1e3, 1e2);

	setup.run();

}

TEST_F(DCDGTest, testTDynamicLabelPropagationStrategyIsolate) {

	DynamicCommunityDetector* dynPLP = new TDynamicLabelPropagation<Isolate>();
	INFO("created algorithm: " << dynPLP->toString());

	DynamicGraphSource* dynGen = new DynamicBarabasiAlbertGenerator(1);

	std::vector<DynamicCommunityDetector*> detectors = {dynPLP};
	DynCDSetup setup(*dynGen, detectors, 1e3, 1e2);

	setup.run();

	Graph G = setup.getGraphCopy();

	for (std::vector<Clustering> clusteringSequence : setup.results) {
		Clustering last = clusteringSequence.back();
		EXPECT_TRUE(last.isProper(G)) << "final clustering in the sequence should be a proper clustering of G";
	}


}

TEST_F(DCDGTest, testTDynamicLabelPropagationStrategyIsolateNeighbors) {

	DynamicCommunityDetector* dynPLP = new TDynamicLabelPropagation<IsolateNeighbors>();
	INFO("created algorithm: " << dynPLP->toString());

	DynamicGraphSource* dynGen = new DynamicBarabasiAlbertGenerator(1);

	std::vector<DynamicCommunityDetector*> detectors = {dynPLP};
	DynCDSetup setup(*dynGen, detectors, 1e3, 1e2);

	setup.run();

	Graph G = setup.getGraphCopy();

	for (std::vector<Clustering> clusteringSequence : setup.results) {
		Clustering last = clusteringSequence.back();
		EXPECT_TRUE(last.isProper(G)) << "final clustering in the sequence should be a proper clustering of G";
	}

}

TEST_F(DCDGTest, testDynamicEnsembleWithTDynamicLabelPropagation) {
	DynamicGraphSource* dynGen = new DynamicBarabasiAlbertGenerator(1);
	DynamicCommunityDetector* dynLP1 = new TDynamicLabelPropagation<Isolate>();
	DynamicCommunityDetector* dynLP2 = new TDynamicLabelPropagation<IsolateNeighbors>();
	Clusterer* PLM = new Louvain();

	DynamicEnsemble* ensemble = new DynamicEnsemble();
	ensemble->addBaseAlgorithm(*dynLP1);
	ensemble->addBaseAlgorithm(*dynLP2);
	ensemble->setFinalAlgorithm(*PLM);

	HashingOverlapper overlapAlgo;
	ensemble->setOverlapper(overlapAlgo);

	std::vector<DynamicCommunityDetector*> detectors = { ensemble };
	DynCDSetup setup(*dynGen, detectors, 1e2, 10);

	setup.run();

	Graph G = setup.getGraphCopy();
	for (std::vector<Clustering> clusteringSequence : setup.results) {
		Clustering last = clusteringSequence.back();
		EXPECT_TRUE(last.isProper(G)) << "final clustering in the sequence should be a proper clustering of G";
	}


}


TEST_F(DCDGTest, testSetupWithStatic) {
	DynamicGraphSource* dynGen = new DynamicBarabasiAlbertGenerator(1);
	DynamicCommunityDetector* dynLP1 = new TDynamicLabelPropagation<Isolate>();


	std::vector<DynamicCommunityDetector*> detectors = { dynLP1 };
	DynCDSetup setup(*dynGen, detectors, 1e2, 10);

	Clusterer* staticAlgo = new LabelPropagation();
	setup.setStatic(staticAlgo);

	setup.run();

	Graph G = setup.getGraphCopy();
	for (std::vector<Clustering> clusteringSequence : setup.results) {
		Clustering last = clusteringSequence.back();
		EXPECT_TRUE(last.isProper(G)) << "final clustering in the sequence should be a proper clustering of G";
	}


}

TEST_F(DCDGTest, tryArxivEval) {

	DynamicCommunityDetector* dynPLP = new TDynamicLabelPropagation<Isolate>();
	INFO("created algorithm: " << dynPLP->toString());

	DynamicDGSParser* dynGen = new DynamicDGSParser("/Users/forigem/KIT/arXivSpider/src/arxivspider/conference-dataset/conference-cs-paper.dgs");

	std::vector<DynamicCommunityDetector*> detectors = {dynPLP};
	DynCDSetup setup(*dynGen, detectors, 1e10, 1e2);

	setup.run();

	Graph G = setup.getGraphCopy();
	INFO("Still alive before the loop");
	Clustering last;
	for (std::vector<Clustering> clusteringSequence : setup.results) {
		last = clusteringSequence.back();
	}
	INFO("Still alive before the clusterings");
	dynGen -> evaluateClusterings("clustering-output.txt", last);

//	std::string path = "output-clusts.txg";
	INFO("Is proper: " << last.isProper(G));
	INFO("Number of nodes: " << G.numberOfNodes());
	INFO("Number of edges: " << G.numberOfEdges());

}

} /* namespace NetworKit */

#endif /*NOGTEST */
