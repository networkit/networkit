/*
 * GeneratorsTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "GeneratorsGTest.h"

namespace NetworKit {

GeneratorsGTest::GeneratorsGTest() {
	// TODO Auto-generated constructor stub

}

GeneratorsGTest::~GeneratorsGTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(GeneratorsGTest, testDynamicBarabasiAlbertGenerator) {

	Graph G(0); // empty graph
	GraphEventProxy proxy(G);

	DynamicGraphGenerator* gen = new DynamicBarabasiAlbertGenerator(proxy, 2);

	gen->initializeGraph();

	EXPECT_EQ(2, G.numberOfNodes()) << "initially the generator creates two connected nodes";
	EXPECT_EQ(1, G.numberOfEdges()) << "initially the generator creates two connected nodes";

	count n = 100;

	gen->generateWhile([&]() {
				return ( G.numberOfNodes() < n );
			});

	EXPECT_EQ(n, G.numberOfNodes());
	DEBUG("m = " << G.numberOfEdges());

	// resume generator
	gen->generateWhile([&]() {
		return (G.numberOfNodes() < 2 * n);
	});
	EXPECT_EQ(2 * n, G.numberOfNodes());
}


TEST_F(GeneratorsGTest, viewDynamicBarabasiAlbertGenerator) {
	Graph G(0); // empty graph
	GraphEventProxy proxy(G);
	DynamicGraphGenerator* gen = new DynamicBarabasiAlbertGenerator(proxy, 2);
	gen->initializeGraph();
	count n = 42;
	gen->generateWhile([&]() {
				return ( G.numberOfNodes() < n );
			});
	METISGraphWriter writer;
	writer.write(G, "output/BATest.graph");
}


TEST_F(GeneratorsGTest, testStaticPubWebGenerator) {
	count n = 800;
	count numCluster = 40;
	count maxNumNeighbors = 20;
	float rad = 0.125;

	PubWebGenerator gen(n, numCluster, rad, maxNumNeighbors);
	Graph G = gen.generate();

	EXPECT_EQ(n, G.numberOfNodes()) << "number of generated nodes";

	LabelPropagation lp;
	Clustering clustering = lp.run(G);

	// output to EPS file
	PostscriptWriter psWriter(G, true);
	psWriter.write(clustering, "output/pubweb-lp-cluster.eps");
}


TEST_F(GeneratorsGTest, tryBTERGenerator) {
	std::vector<count> degreeDistribution { 0, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	std::vector<double> clusteringCoefficients {0.0, 0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	DEBUG("construct BTERGenerator");
	BTERGenerator bter(degreeDistribution, clusteringCoefficients, 1.0);
	std::pair<count, count> nm = bter.desiredGraphSize();
	DEBUG("call BTERGenerator");
	Graph G = bter.generate();

	EXPECT_EQ(nm.first, G.numberOfNodes());
	EXPECT_EQ(nm.second, G.numberOfEdges());

	METISGraphWriter writer;
	writer.write(G, "output/BTERTest.graph");
}


TEST_F(GeneratorsGTest, tryBTERGeneratorWithPowerLawDistribution) {
	std::vector<double> clusteringCoefficients {0.0, 0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	std::vector<count> degreeDistribution = BTERGenerator::generatePowerLawDegreeDistribution(14, 0.5);
	DEBUG("construct BTERGenerator");
	BTERGenerator bter(degreeDistribution, clusteringCoefficients, 1.0);
	std::pair<count, count> nm = bter.desiredGraphSize();
	DEBUG("call BTERGenerator");
	Graph G = bter.generate();

	EXPECT_EQ(nm.first, G.numberOfNodes());
	EXPECT_EQ(nm.second, G.numberOfEdges());

	METISGraphWriter writer;
	writer.write(G, "output/BTERTest.graph");
}

TEST_F(GeneratorsGTest, testStaticBarabasiAlbertGenerator) {
	count k = 3;
	count nMax = 100;
	count n0 = 3;

	StaticBarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());

	G = BarabasiAlbert.generate();
	EXPECT_FALSE(G.isEmpty());

	EXPECT_EQ(nMax, G.numberOfNodes());
	EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());



}

TEST_F(GeneratorsGTest, generatetStaticBarabasiAlbertGeneratorGraph) {
		count k = 3;
		count nMax = 1000;
		count n0 = 3;

		StaticBarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);

		Graph G = BarabasiAlbert.generate();
		GraphIO io;
		io.writeAdjacencyList(G, "output/"
				"BarabasiGraph.txt");
}


} /* namespace NetworKit */

