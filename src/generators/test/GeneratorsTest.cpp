/*
 * GeneratorsTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#include "GeneratorsTest.h"

namespace NetworKit {

GeneratorsTest::GeneratorsTest() {
	// TODO Auto-generated constructor stub

}

GeneratorsTest::~GeneratorsTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(GeneratorsTest, testDynamicBarabasiAlbertGenerator) {

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

TEST_F(GeneratorsTest, testStaticPubWebGenerator) {
	count n = 800;
	count numCluster = 40;
	count maxNumNeighbors = 20;
	float rad = 0.125;

	PubWebGenerator gen(n, numCluster, rad, maxNumNeighbors);
	Graph G = gen.generate();

	EXPECT_EQ(n, G.numberOfNodes()) << "number of generated nodes";

	ClusteringGenerator clusterGen;
	Clustering clustering = clusterGen.makeRandomClustering(G, numCluster);

	// output to EPS file
	PostscriptWriter psWriter(G, true);
	psWriter.write(clustering, "pubweb-random-cluster.eps");
}

} /* namespace NetworKit */

