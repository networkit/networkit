/*
Dy * GeneratorsTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef NOGTEST

#include "GeneratorsGTest.h"

namespace NetworKit {

GeneratorsGTest::GeneratorsGTest() {
	// TODO Auto-generated constructor stub

}

GeneratorsGTest::~GeneratorsGTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(GeneratorsGTest, testDynamicBarabasiAlbertGeneratorSingleStep) {
	count k = 2; // number of edges added per node
	DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(k);
	GraphEventProxy* Gproxy = gen->newGraph();
	Graph* G = Gproxy->G;

	gen->initializeGraph();

	count nPre = G->numberOfNodes();
	count mPre = G->numberOfEdges();
	EXPECT_EQ(k, nPre) << "graph should have been initialized to k nodes";
	EXPECT_EQ(k - 1, mPre) << "graph should have been initialized to a path of k nodes which means k-1 edges";

	// perform single preferential attachment step
	gen->generate();

	count nPost = G->numberOfNodes();
	count mPost = G->numberOfEdges();
	EXPECT_EQ(nPre + 1, nPost) << "one more node should have been added";
	EXPECT_EQ(mPre + k, mPost) << "k edges should have been added";

	delete G;
}

TEST_F(GeneratorsGTest, testDynamicBarabasiAlbertGenerator) {


	DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(2);

	GraphEventProxy* Gproxy = gen->newGraph();
	Graph* G = Gproxy->G;

	gen->initializeGraph();

	EXPECT_EQ(2, G->numberOfNodes()) << "initially the generator creates two connected nodes";
	EXPECT_EQ(1, G->numberOfEdges()) << "initially the generator creates two connected nodes";

	count n = 100;

	gen->generateWhile([&]() {
				return ( G->numberOfNodes() < n );
			});

	EXPECT_EQ(n, G->numberOfNodes());
	DEBUG("m = " << G->numberOfEdges());

	// resume generator

	gen->generateWhile([&]() {
		return (G->numberOfNodes() < 2 * n);
	});
	EXPECT_EQ(2 * n, G->numberOfNodes());
}


TEST_F(GeneratorsGTest, viewDynamicBarabasiAlbertGenerator) {
	DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(2);
	GraphEventProxy* Gproxy = gen->newGraph();
	Graph* G = Gproxy->G;
	gen->initializeGraph();
	count n = 42;
	gen->generateWhile([&]() {
				return ( G->numberOfNodes() < n );
			});
	METISGraphWriter writer;
	writer.write(*G, "output/BATest.graph");

	delete G;
}


TEST_F(GeneratorsGTest, testStaticPubWebGenerator) {
	count n = 1000;
	count numCluster = 10;
	count maxNumNeighbors = 8;
	float rad = 0.15;

	PubWebGenerator gen(n, numCluster, rad, maxNumNeighbors);
	Graph G = gen.generate();

	EXPECT_EQ(n, G.numberOfNodes()) << "number of generated nodes";

	// check degree
	G.forNodes([&](node v) {
		EXPECT_LE(G.degree(v), maxNumNeighbors) << "maximum degree";
	});

	// 1-clustering
	ClusteringGenerator clusterGen;
	Clustering oneClustering = clusterGen.makeOneClustering(G);

	// output to EPS file
	PostscriptWriter psWriter(G, true);
	psWriter.write(oneClustering, "output/pubweb.eps");

	// clustering
	PLM2 clusterAlgo;
	Clustering clustering = clusterAlgo.run(G);
	PostscriptWriter psWriter2(G, true);
	psWriter2.write(clustering, "output/pubweb-clustered-plm2.eps");

	Modularity mod;
	double modVal = mod.getQuality(clustering, G);
	EXPECT_GE(modVal, 0.3) << "modularity of clustering";
}


// FIXME: segmentation fault
TEST_F(GeneratorsGTest, tryDynamicPubWebGenerator) {

	count numInitialNodes = 300;
	count numberOfDenseAreas = 10;
	float neighborhoodRadius = 0.125;
	count maxNumberOfNeighbors = 16;
	count numIterations = 10;


	DynamicGraphSource* gen = new DynamicPubWebGenerator(numInitialNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors);

	GraphEventProxy* Gproxy = gen->newGraph();
	Graph* G = Gproxy->G;

	TRACE("before init graph");
	gen->initializeGraph();
	TRACE("after init graph");

	EXPECT_EQ(numInitialNodes, G->numberOfNodes()) << "initial number of nodes";

	TRACE("m = " << G->numberOfEdges());

	for (index iter = 0; iter < numIterations; ++iter) {
		gen->generate();
		TRACE("m = " << G->numberOfEdges());

		PostscriptWriter psWriter(*G, true);
		char filename[20];
		assert(iter < 10);
		sprintf(filename, "output/pubweb-%i.eps", int(iter));
		psWriter.write(filename);
	}
}



TEST_F(GeneratorsGTest, testBarabasiAlbertGenerator) {
	count k = 3;
	count nMax = 100;
	count n0 = 3;

	BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());

	G = BarabasiAlbert.generate();
	EXPECT_FALSE(G.isEmpty());

	EXPECT_EQ(nMax, G.numberOfNodes());
	EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());



}

TEST_F(GeneratorsGTest, generatetBarabasiAlbertGeneratorGraph) {
		count k = 3;
		count nMax = 1000;
		count n0 = 3;

		BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);

		Graph G = BarabasiAlbert.generate();
		GraphIO io;
		io.writeAdjacencyList(G, "output/"
				"BarabasiGraph.txt");
}


} /* namespace NetworKit */

#endif /*NOGTEST */

