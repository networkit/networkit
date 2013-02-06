/*
 * GraphGeneratorGTest.cpp
 *
 *  Created on: 31.01.2013
 *      Author: cls
 */

#include "GraphGeneratorGTest.h"

namespace EnsembleClustering {

GraphGeneratorGTest::GraphGeneratorGTest() {
	// TODO Auto-generated constructor stub

}

GraphGeneratorGTest::~GraphGeneratorGTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(GraphGeneratorGTest, testBarabasiAlbert) {
	int64_t n = 5000;
	int64_t k = 2;
	GraphGenerator graphGen;
	// Graph G = graphGen.makePreferentialAttachmentGraph(n, k);

	// INFO("m = " << G.numberOfEdges());

	GraphIO graphio;
	INFO("writing graph to file");
	// graphio.writeAdjacencyList(G, "sandbox/testBarabasiAlbert.graph");

	EXPECT_TRUE(false) << "TODO: fix preferential attachment";
}


} /* namespace EnsembleClustering */
