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

	DynamicGraphGenerator* gen = new DynamicBarabasiAlbertGenerator(proxy);

	gen->initializeGraph();

	EXPECT_EQ(2, G.numberOfNodes()) << "initially the generator creates two connected nodes";

	count n = 100;

	gen->generate([&]() {
				return ( G.numberOfNodes() == n );
			});

	EXPECT_EQ(n, G.numberOfNodes());

	DEBUG("m = " << G.numberOfEdges());
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


TEST_F(GeneratorsGTest, testBTERGenerator) {
	std::vector<count> degreeDistribution { 0, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	std::vector<double> clusteringCoefficients {0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	DEBUG("construct BTERGenerator");
	BTERGenerator bter(degreeDistribution, clusteringCoefficients, 1.0);
	DEBUG("call BTERGenerator");
	bter.generate();
}

} /* namespace NetworKit */

