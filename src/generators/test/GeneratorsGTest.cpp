/*
Dy * GeneratorsTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef NOGTEST

#include "GeneratorsGTest.h"

#include "../DynamicPathGenerator.h"
#include "../ForestFireGenerator.h"

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
	DEBUG("m = " , G->numberOfEdges());

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
	count n = 1800;
	count numCluster = 24;
	count maxNumNeighbors = 36;
	float rad = 0.075;

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
	psWriter.write(clustering, "output/pubweb-clustered-plm2.eps");

	Modularity mod;
	double modVal = mod.getQuality(clustering, G);
	EXPECT_GE(modVal, 0.2) << "modularity of clustering";
	DEBUG("Modularity of clustering: " , modVal);
	DEBUG("Total edge weight: " , G.totalEdgeWeight());
}


TEST_F(GeneratorsGTest, testDynamicPubWebGenerator) {
	count nSteps = 10;
	count n = 400;
	count numCluster = 12;
	count maxNumNeighbors = 36;
	float rad = 0.125;

	DynamicPubWebGenerator dynGen(n, numCluster, rad, maxNumNeighbors, false);
	Graph G = dynGen.getGraph();
	GraphUpdater gu(G);
	std::vector<GraphEvent> stream;

	// static clustering algorithm for better visual output
	PostscriptWriter psWriter(G, true);
	psWriter.write("output/pubweb-0000.eps");

	for (index i = 1; i <= nSteps; ++i) {
		stream = dynGen.generate(1);
		DEBUG("updating graph");
		gu.update(stream);
		G.initCoordinates();

		DEBUG("updated graph, new (n, m) = (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");

		// update coordinates
		std::map<node, Point<float> > newCoordinates = dynGen.getNewCoordinates();
		for (std::map<node, Point<float> >::iterator iter = newCoordinates.begin();
				iter != newCoordinates.end(); ++iter) {
			node v = iter->first;
			Point<float> p = iter->second;
			G.setCoordinate(v, p);
		}

		// output in clustered way for better visual inspection
		char path[23];
		sprintf(path, "output/pubweb-%04llu.eps", i);
		TRACE("path: " , path);
		PostscriptWriter psWriter(G, true);
		psWriter.write(path);
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


TEST_F(GeneratorsGTest, testDynamicPathGenerator) {
	DynamicPathGenerator gen;
	auto stream = gen.generate(42);
	for (auto ev : stream) {
		TRACE(ev.toString());
	}
}


TEST_F(GeneratorsGTest, testErdosRenyiGenerator) {
	count n = 2000;
	double p = 1.5 * (log(n) / (double) n);

	ErdosRenyiGenerator generator(n, p);
	Graph G = generator.generate();
	EXPECT_EQ(n, G.numberOfNodes());
	EXPECT_FALSE(G.isEmpty());

	count nPairs = (n * (n-1)) / 2;
	count nEdges = G.numberOfEdges();
	EXPECT_GE(nEdges, 0.75 * p * nPairs);
	EXPECT_LE(nEdges, 1.25 * p * nPairs);

	DEBUG("Number of edges with probability " , p , " (actual/expected): " , nEdges , " / " , (nPairs * p));
}

TEST_F(GeneratorsGTest, testChungLuGenerator) {
	count n = 200;
	count maxDegree = n / 8;
	std::vector<count> sequence(n);
	count expVolume = 0;
	count actualVolume = 0;

	// fill sequence with random values (this is not power-law, of course!)
	for (index i = 0; i < n; ++i) {
		sequence[i] = rand() % maxDegree;
		expVolume += sequence[i];
	}

	ChungLuGenerator gen(sequence);
	Graph G = gen.generate();

	EXPECT_EQ(n, G.numberOfNodes());
	G.forNodes([&](node v) {
		actualVolume += G.degree(v);
	});

	INFO("expected volume: ", expVolume, ", actual volume: ", actualVolume);
}


TEST_F(GeneratorsGTest, testForestFireGenerator) {
	ForestFireGenerator ffg(0.5);
	ffg.generate(10);
}

} /* namespace NetworKit */

#endif /*NOGTEST */

